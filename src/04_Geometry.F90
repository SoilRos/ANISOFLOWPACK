MODULE ANISOFLOW_Geometry

! ANISOFLOW_Geometry it's a module that contains routines to manage geometry variables.

    IMPLICIT NONE

CONTAINS

 !  - GetGeometry: It's a routine that fills a Geometry data structure with input files provided by
 !                 the user.
 !    > OUT: Gmtry, ierr.
 !      + Gmtry: It's a Geometry data structure filled with input files provided by the user.
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.
 !    > NOTES: The Geometry is filled in three stages: to create a Data Manager (DataMngr) to control
 !             information related to a regular rectangular grid, a grid indformation and, a Topology 
 !             that describes the geometry with identifiers.

SUBROUTINE GetGeometry(Comm,Gmtry,ierr)

    USE ANISOFLOW_Interface,    ONLY : GetVerbose
    USE ANISOFLOW_Types,        ONLY : Geometry

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    MPI_Comm,INTENT(IN)             :: Comm
    TYPE(Geometry),INTENT(OUT)      :: Gmtry

    PetscBool                       :: Verbose
    CHARACTER(LEN=200)              :: EventName,ClassName
    PetscLogEvent                   :: Event
    PetscClassId                    :: ClassID
    PetscLogDouble                  :: EventFlops

    ClassName="Geometry"
    CALL PetscClassIdRegister(ClassName,ClassID,ierr)
    EventName="GetGeometry"
    CALL PetscLogEventRegister(EventName,ClassID,Event,ierr)
    CALL PetscLogEventBegin(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,"["//ADJUSTL(TRIM(EventName))//" Event] Inizialited\n",ierr)
    
    CALL GetDataMngr(Comm,Gmtry%DataMngr,Gmtry%Scale,ierr)
    CALL GetGrid(Comm,Gmtry%DataMngr,Gmtry%x,Gmtry%y,Gmtry%z,Gmtry%Scale,ierr)
    CALL GetTopology(Comm,Gmtry%DataMngr,Gmtry%Tplgy,Gmtry%DirichIS,Gmtry%SourceIS,Gmtry%CauchyIS,Gmtry%Scale,ierr)
    
    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,"["//ADJUSTL(TRIM(EventName))//" Event] Finalized\n",ierr)
    
    CALL PetscLogFlops(EventFlops,ierr)
    CALL PetscLogEventEnd(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

END SUBROUTINE GetGeometry

 !  - GetDataMngr: It's a routine that creates and fills the information related with a regular 
 !                 rectangular grid.
 !    > OUT: DataMngr, ierr.
 !      + DataMngr: It's a DMDA PETSc structure that stores the information related with a regular 
 !                 rectangular grid.
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.
 !    > NOTES: The DataMngr takes the size of the domain and assign a subdomain to each processor.

SUBROUTINE GetDataMngr(Comm,DataMngr,Scale,ierr)

    USE ANISOFLOW_Types,        ONLY : InputTypeVar
    USE ANISOFLOW_Interface,    ONLY : GetInputType,GetVerbose

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscdm.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    MPI_Comm,INTENT(IN)             :: Comm
    PetscInt,INTENT(IN)             :: Scale
    DM,INTENT(OUT)                  :: DataMngr

    TYPE(InputTypeVar)              :: InputType

    IF (Scale.EQ.1) THEN
        CALL GetInputType(InputType,ierr)
    END IF

    ! InputType define the type of file that is provided.
    !   1: Defined by Blessent. An example is provided in "../ex/Blessent/in/tsim_USMH.asc"
    !   2: Defined by Perez. An example is provided in "../ex/Perez/in/sanpck.domnRST"

    IF (InputType%Gmtry.EQ.1) THEN
        CALL GetDataMngr_1(Comm,DataMngr,Scale,ierr)
    ELSE IF (InputType%Gmtry.EQ.2) THEN
        CALL GetDataMngr_2(Comm,DataMngr,Scale,ierr)
    ELSE
        CALL PetscSynchronizedPrintf(Comm, &
            & "[ERROR] Geometry InputType wrong\n",ierr)
        STOP
    END IF

END SUBROUTINE GetDataMngr

 !  - GetDataMngr_1: It's a routine that creates and fills the information related with a regular 
 !                   rectangular grid when InputType%Gmtry=1.
 !    > OUT: DataMngr, ierr.
 !      + DataMngr: It's a DMDA PETSc structure that stores the information related with a regular 
 !                 rectangular grid.
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.
 !    > NOTES: The DataMngr takes the size of the domain and assign a subdomain to each processor.

SUBROUTINE GetDataMngr_1(Comm,DataMngr,Scale,ierr)

    USE ANISOFLOW_Types,        ONLY : RunOptionsVar
    USE ANISOFLOW_Interface,    ONLY : GetInputDir,GetInputFileGmtry,GetRunOptions,GetVerbose

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    MPI_Comm,INTENT(IN)             :: Comm
    PetscInt,INTENT(IN)             :: Scale
    DM,INTENT(OUT)                  :: DataMngr

    PetscMPIInt                     :: process
    PetscInt                        :: u,widthG(3)
    CHARACTER(LEN=200)              :: InputDir,InputFileGmtry,Route
    TYPE(RunOptionsVar)             :: RunOptions
    DMDAStencilType                 :: Stencil
    PetscBool                       :: Verbose

    PARAMETER(u=01)

    CALL GetVerbose(Verbose,ierr)

    ! It obtains the route to open a geometry file.
    IF (Scale.EQ.1) THEN
        CALL GetInputDir(InputDir,ierr)
        CALL GetInputFileGmtry(InputFileGmtry,ierr)
    END IF

    ! It obtains run options.
    CALL GetRunOptions(RunOptions,ierr)

    CALL MPI_Comm_rank(Comm,process,ierr)

    ! It obtains the global size of the domain on the first processor.
    IF (process.EQ.0) THEN
        Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFileGmtry))
        OPEN(u,FILE=TRIM(Route),STATUS='OLD',ACTION='READ')

        ! It gets the size of the domain 
        READ(u, '((I10),(I10),(I10))')widthG(1),widthG(2),widthG(3)
        CLOSE(u)
    END IF

    ! It broadcasts the global size to other processors.
    CALL MPI_Bcast(widthG,3,MPI_INT,0, Comm,ierr)

    ! It decides the stencil shape depending on the scheme used.
    IF (RunOptions%Scheme.EQ.1) THEN
        Stencil=DMDA_STENCIL_STAR
    ELSEIF (RunOptions%Scheme.EQ.2) THEN
        Stencil=DMDA_STENCIL_BOX
    ELSE
        CALL PetscSynchronizedPrintf(Comm,                         &
            & "[ERROR] Run_options_scheme command must be an integer between 1 and 2\n",ierr)
        STOP
    END IF

    ! It creates the Data Manager to Distributed Arrays by the information provided.
    IF (widthG(3).NE.1) THEN
        CALL DMDACreate3d(Comm,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,      &
            & DM_BOUNDARY_GHOSTED,Stencil,widthG(1),widthG(2),widthG(3),     &
            & PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,PETSC_NULL_INTEGER,       &
            & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,DataMngr,ierr)
    ELSE
        CALL DMDACreate2d(Comm,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,      &
            & Stencil,widthG(1),widthG(2),PETSC_DECIDE,PETSC_DECIDE,1,1,PETSC_NULL_INTEGER,       &
            & PETSC_NULL_INTEGER,DataMngr,ierr)
    END IF

    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,"[GetGeometry Event] Data Manager was satisfactorily created\n",ierr)

END SUBROUTINE GetDataMngr_1

SUBROUTINE GetDataMngr_2(Comm,DataMngr,Scale,ierr)

    USE ANISOFLOW_Types,        ONLY : RunOptionsVar
    USE ANISOFLOW_Interface,    ONLY : GetInputDir,GetInputFileGmtry,GetRunOptions,GetVerbose

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    MPI_Comm,INTENT(IN)             :: Comm
    PetscInt,INTENT(IN)             :: Scale
    DM,INTENT(OUT)                  :: DataMngr

    PetscMPIInt                     :: process
    PetscInt                        :: u,widthG(3)
    CHARACTER(LEN=200)              :: InputDir,InputFileGmtry,Route,aux
    TYPE(RunOptionsVar)             :: RunOptions
    DMDAStencilType                 :: Stencil
    PetscBool                       :: Verbose

    PARAMETER(u=01)

    CALL GetVerbose(Verbose,ierr)

    ! It obtains the route to open a geometry file.
    IF (Scale.EQ.1) THEN
        CALL GetInputDir(InputDir,ierr)
        CALL GetInputFileGmtry(InputFileGmtry,ierr)
    END IF

    ! It obtains run options.
    CALL GetRunOptions(RunOptions,ierr)

    ! It tag every processor from 0 to n-1
    CALL MPI_Comm_rank(Comm,process,ierr)

    ! It obtains the global size of the domain on the first processor.
    IF (process.EQ.0) THEN
        Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFileGmtry))
        OPEN(u,FILE=TRIM(Route),STATUS='OLD',ACTION='READ')

        ! It gets the size of the domain 
        READ(u, '((A10),(I10))')aux,widthG(1)
        READ(u, '((A10),(I10))')aux,widthG(2)
        READ(u, '((A10),(I10))')aux,widthG(3)
        CLOSE(u)
    END IF

    ! It broadcasts the global size to other processors.
    CALL MPI_Bcast(widthG,3,MPI_INT, 0, Comm,ierr)

    ! It decides the stencil shape depending on the scheme used.
    IF (RunOptions%Scheme.EQ.1) THEN
        Stencil=DMDA_STENCIL_STAR
    ELSEIF (RunOptions%Scheme.EQ.2) THEN
        Stencil=DMDA_STENCIL_BOX
    ELSE
        CALL PetscSynchronizedPrintf(Comm,                         &
            & "[ERROR] Run_options_scheme command must be an integer between 1 and 2\n",ierr)
        STOP
    END IF

    ! It creates the Data Manager to Distributed Arrays by the information provided.
    IF (widthG(3).NE.1) THEN
        CALL DMDACreate3d(Comm,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,      &
            & DM_BOUNDARY_GHOSTED,Stencil,widthG(1),widthG(2),widthG(3),     &
            & PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,PETSC_NULL_INTEGER,       &
            & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,DataMngr,ierr)
    ELSE
        CALL DMDACreate2d(Comm,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,      &
            & Stencil,widthG(1),widthG(2),PETSC_DECIDE,PETSC_DECIDE,1,1,PETSC_NULL_INTEGER,       &
            & PETSC_NULL_INTEGER,DataMngr,ierr)
    END IF

    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,"[GetGeometry Event] Data Manager was satisfactorily created\n",ierr)

END SUBROUTINE GetDataMngr_2

 !  - GetGrid: It's a routine that creates and fills the information related with a a coordinates of
 !             a regular rectangular grid.
 !    > IN: DataMngr
 !      + DataMngr: It's a DMDA PETSc structure that has stored the information related with a 
 !                 regular rectangular grid.
 !    > OUT: x, y, z, ierr.
 !      + x,y,z: It's an Array that has de coordinates of the grid on each direction. It's stored on 
 !               every processor 
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.

SUBROUTINE GetGrid(Comm,DataMngr,x,y,z,Scale,ierr)

    USE ANISOFLOW_Types,        ONLY : InputTypeVar
    USE ANISOFLOW_Interface,    ONLY : GetInputType,GetVerbose

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    MPI_Comm,INTENT(IN)                 :: Comm
    PetscInt,INTENT(IN)                 :: Scale
    DM,INTENT(IN)                       :: DataMngr
    Vec,INTENT(OUT)                     :: x,y,z

    TYPE(InputTypeVar)                  :: InputType

    IF (Scale.EQ.1) THEN
        CALL GetInputType(InputType,ierr)
    END IF

    ! InputType define the type of file that is provided.
    !   1: Defined by default. Default grid used has DX=DY=DZ=1.0
    IF (InputType%Gmtry.EQ.1) THEN
        CALL GetGrid_1(Comm,DataMngr,x,y,z,Scale,ierr)
    ELSE IF (InputType%Gmtry.EQ.2) THEN
        CALL GetGrid_2(Comm,DataMngr,x,y,z,Scale,ierr)
    ELSE        
        CALL PetscSynchronizedPrintf(Comm, &
            & "[ERROR] Geometry InputType wrong\n",ierr)
        STOP
    END IF

END SUBROUTINE GetGrid

 !  - GetGrid_1: It's a routine that creates and fills the information related with a a coordinates of
 !               a regular rectangular grid when InputType%Gmtry=1
 !    > IN: DataMngr
 !      + DataMngr: It's a DMDA PETSc structure that has stored the information related with a 
 !                 regular rectangular grid.
 !    > OUT: x, y, z, ierr.
 !      + x,y,z: It's an Array that has de coordinates of the grid on each direction. It's stored on 
 !               every processor 
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.

SUBROUTINE GetGrid_1(Comm,DataMngr,x,y,z,Scale,ierr)

    USE ANISOFLOW_Interface, ONLY : GetInputDir,GetInputFileGmtry,GetVerbose

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    MPI_Comm,INTENT(IN)                 :: Comm
    PetscInt,INTENT(IN)                 :: Scale
    DM,INTENT(IN)                       :: DataMngr
    Vec,INTENT(OUT)                     :: x,y,z

    CHARACTER(LEN=200)                  :: InputDir,InputFileGmtry
    PetscInt                            :: widthG(3),size,i
    PetscReal                           :: Value
    PetscBool                           :: Verbose

    CALL GetVerbose(Verbose,ierr)

    ! It obtains the route to open a geometry file.
    IF (Scale.EQ.1) THEN
        CALL GetInputDir(InputDir,ierr)
        CALL GetInputFileGmtry(InputFileGmtry,ierr)
    END IF

    ! It gets the global size from the geometry data manager.
    CALL DMDAGetInfo(DataMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
        & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,ierr)

    ! It create the vectors
    size=widthG(1)+1
    CALL VecCreateSeq(PETSC_COMM_SELF,size,x,ierr)

    size=widthG(2)+1
    CALL VecCreateSeq(PETSC_COMM_SELF,size,y,ierr)

    size=widthG(3)+1
    CALL VecCreateSeq(PETSC_COMM_SELF,size,z,ierr)

    
    DO i=0,widthG(1)
        Value=REAL(i-1)
        CALL VecSetValue(x,i,Value,INSERT_VALUES,ierr)
    END DO
    DO i=0,widthG(2)
        Value=REAL(i-1)
        CALL VecSetValue(y,i,Value,INSERT_VALUES,ierr)
    END DO
    DO i=0,widthG(3)
        Value=REAL(i-1)
        CALL VecSetValue(z,i,Value,INSERT_VALUES,ierr)
    END DO
    ! Default DX=DY=DZ=1.0
    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,"[GetGeometry Event] WARNING: Grid System wasn't provided. Default grid used has DX=DY=DZ=1.0\n",ierr)

    CALL VecAssemblyBegin(x,ierr)
    CALL VecAssemblyEnd(x,ierr)

    CALL VecAssemblyBegin(y,ierr)
    CALL VecAssemblyEnd(y,ierr)

    CALL VecAssemblyBegin(z,ierr)
    CALL VecAssemblyEnd(z,ierr)

    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,"[GetGeometry Event] Grid System was satisfactorily created\n",ierr)

END SUBROUTINE GetGrid_1

SUBROUTINE GetGrid_2(Comm,DataMngr,x,y,z,Scale,ierr)
    
    USE ANISOFLOW_Interface, ONLY : GetInputDir,GetInputFileGmtry,GetVerbose

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    MPI_Comm,INTENT(IN)                 :: Comm
    PetscInt,INTENT(IN)                 :: Scale
    DM,INTENT(IN)                       :: DataMngr
    Vec,INTENT(OUT)                     :: x,y,z

    PetscMPIInt                         :: process
    CHARACTER(LEN=200)                  :: InputDir,InputFileGmtry,Route
    PetscInt                            :: widthG(3),size,i,u
    PetscReal                           :: Value
    PetscBool                           :: Verbose

    PARAMETER(u=01)

    CALL GetVerbose(Verbose,ierr)

    ! It obtains the route to open a geometry file.
    IF (Scale.EQ.1) THEN
        CALL GetInputDir(InputDir,ierr)
        CALL GetInputFileGmtry(InputFileGmtry,ierr)
    END IF


    ! It gets the global size from the geometry data manager.
    CALL DMDAGetInfo(DataMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
        & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,ierr)

    ! It create the vectors
    size=widthG(1)+1
    CALL VecCreateSeq(PETSC_COMM_SELF,size,x,ierr)

    size=widthG(2)+1
    CALL VecCreateSeq(PETSC_COMM_SELF,size,y,ierr)

    size=widthG(3)+1
    CALL VecCreateSeq(PETSC_COMM_SELF,size,z,ierr)

    ! It tag every processor from 0 to n-1
    CALL MPI_Comm_rank(Comm,process,ierr)

    ! It obtains the global size of the domain on the first processor.
    IF (process.EQ.0) THEN
        Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFileGmtry))
        OPEN(u,FILE=TRIM(Route),STATUS='OLD',ACTION='READ')
        READ(u,*)
        READ(u,*)
        READ(u,*)
    END IF

    DO i=0,widthG(1)
        IF (process.EQ.0) THEN
            READ(u, '(ES17.11)')Value
        END IF
        ! It broadcasts the global size to other processors.
        CALL MPI_Bcast(Value,1,MPI_DOUBLE, 0, Comm,ierr)
        CALL VecSetValue(x,i,Value,INSERT_VALUES,ierr)
    END DO
    
    DO i=0,widthG(2)
        IF (process.EQ.0) THEN
            READ(u, '(ES17.11)')Value
        END IF
        ! It broadcasts the global size to other processors.
        CALL MPI_Bcast(Value,1,MPI_DOUBLE, 0, Comm,ierr)
        CALL VecSetValue(y,i,Value,INSERT_VALUES,ierr)
    END DO

    DO i=0,widthG(3)
        IF (process.EQ.0) THEN
            READ(u, '(ES17.11)')Value
        END IF
        ! It broadcasts the global size to other processors.
        CALL MPI_Bcast(Value,1,MPI_DOUBLE, 0, Comm,ierr)
        CALL VecSetValue(z,i,Value,INSERT_VALUES,ierr)
    END DO

    IF (process.EQ.0) THEN
        CLOSE(u)
    END IF

    CALL VecAssemblyBegin(x,ierr)
    CALL VecAssemblyEnd(x,ierr)

    CALL VecAssemblyBegin(y,ierr)
    CALL VecAssemblyEnd(y,ierr)

    CALL VecAssemblyBegin(z,ierr)
    CALL VecAssemblyEnd(z,ierr)

    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,"[GetGeometry Event] Grid System was satisfactorily created\n",ierr)

END SUBROUTINE GetGrid_2

 !  - GetTopology: It's a routine that creates and fills the information related to topology.
 !                 It creates a vector and index sets to describe the geometry. 
 !    > IN: DataMngr.
 !      + DataMngr: It's a DMDA PETSc structure that has stored the information related with a 
 !                 regular rectangular grid.
 !    > OUT: Tplgy, DirichIS, NeummanIS, CauchyIS, SourceIS, ierr.
 !      + Tplgy: It's a PETSc vector type produced by DataMngr that contains a topology identifier
 !               in each cell.
 !          0: Inactive cell.
 !          1: Active cell.
 !          2: Dirichlet boundary condition cell.
 !          3: Source Q.
 !          4: Cauchy boundary condition
 !      + DirichIS: It's a PETSc index set that has a map between Dirichlet information and vecs 
 !                  produced by DataMngr.
 !      + SourceIS: It's a PETSc index set that has a map between Source information and vecs 
 !                   roduced by DataMngr.
 !      + CauchyIS: It's a PETSc index set that has a map between Cauchy information and vecs 
 !                  produced by DataMngr.
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.

SUBROUTINE GetTopology(Comm,DataMngr,Tplgy,DirichIS,SourceIS,CauchyIS,Scale,ierr)

    USE ANISOFLOW_Types, ONLY : InputTypeVar
    USE ANISOFLOW_Interface, ONLY : GetInputType,GetVerbose

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    MPI_Comm,INTENT(IN)             :: Comm
    PetscInt,INTENT(IN)             :: Scale
    DM,INTENT(IN)                   :: DataMngr
    Vec,INTENT(OUT)                 :: Tplgy
    IS,INTENT(OUT)                  :: DirichIS,CauchyIS,SourceIS

    TYPE(InputTypeVar)              :: InputType
    PetscInt                        :: SizeDirichIS,SizeSourceIS,SizeCauchyIS

    IF (Scale.EQ.1) THEN
        CALL GetInputType(InputType,ierr)
    END IF

    ! InputType define the type of file that is provided.
    !   1: Default topology, it doesn't need a file. The first border layer is Dirichlet, active in other case.
    IF (InputType%Tplgy.EQ.1) THEN
        CALL GetTopology_1(Comm,DataMngr,Tplgy,DirichIS,SourceIS,CauchyIS,Scale,ierr)
    ELSE IF (InputType%Tplgy.EQ.2) THEN
        CALL GetTopology_2(Comm,DataMngr,Tplgy,DirichIS,SourceIS,CauchyIS,Scale,ierr)
    ELSE IF (InputType%Tplgy.EQ.3) THEN
        CALL GetTopology_3(Comm,DataMngr,Tplgy,DirichIS,SourceIS,CauchyIS,Scale,ierr)
    ELSE
        CALL PetscSynchronizedPrintf(Comm, &
            & "[ERROR] Topology InputType wrong\n",ierr)
        STOP
    END IF

    CALL ISGetLocalSize(DirichIS,SizeDirichIS,ierr)
    CALL ISGetLocalSize(SourceIS,SizeSourceIS,ierr)
    CALL ISGetLocalSize(CauchyIS,SizeCauchyIS,ierr)

    IF ((SizeDirichIS+SizeSourceIS+SizeCauchyIS).EQ.0) THEN
        CALL PetscSynchronizedPrintf(Comm,&
            & "[ERROR] It has to have at least one boundary condition.\n",ierr)
        STOP
    END IF

END SUBROUTINE GetTopology

 !  - GetTopology_1: It's a routine that creates and fills the information related to topology.
 !                   when InputType%Tplgy=1. It creates a vector and index sets to describe 
 !                   the geometry. 
 !    > IN: DataMngr.
 !      + DataMngr: It's a DMDA PETSc structure that has stored the information related with a 
 !                 regular rectangular grid. 
 !    > OUT: Tplgy, DirichIS, NeummanIS, CauchyIS, SourceIS, ierr.
 !      + Tplgy: It's a PETSc vector type produced by DataMngr that contains a topology identifier
 !               in each cell.
 !          0: Inactive cell.
 !          1: Active cell.
 !          2: Dirichlet boundary condition cell.
 !          3: Source Q.
 !          4: Cauchy boundary condition
 !      + DirichIS: It's a PETSc index set that has a map between Dirichlet information and vecs 
 !                  produced by DataMngr.
 !      + SourceIS: It's a PETSc index set that has a map between Source information and vecs 
 !                   roduced by DataMngr.
 !      + CauchyIS: It's a PETSc index set that has a map between Cauchy information and vecs 
 !                  produced by DataMngr.
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.

! from file

SUBROUTINE GetTopology_1(Comm,DataMngr,Tplgy,DirichIS,SourceIS,CauchyIS,Scale,ierr)

    USE ANISOFLOW_Interface, ONLY : GetVerbose,GetInputDir,GetInputFileTplgy
    USE ANISOFLOW_View, ONLY : ViewTopology

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscdmda.h90>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    MPI_Comm,INTENT(IN)             :: Comm
    PetscInt,INTENT(IN)             :: Scale
    DM,INTENT(IN)                   :: DataMngr
    Vec,INTENT(OUT)                 :: Tplgy
    IS,INTENT(OUT)                  :: DirichIS,SourceIS,CauchyIS

    PetscMPIInt                     :: process
    CHARACTER(LEN=200)              :: InputDir,InputFileTplgy,Route,ViewName,EventName
    PetscReal                       :: ValR
    PetscInt                        :: u,i,ValI,widthG(3),BCLenL(3),BCLenG(3)
    Vec                             :: TmpTplgy
    PetscBool                       :: Verbose

    PARAMETER(u=01)

    CALL GetVerbose(Verbose,ierr)

    ! It obtains the route to open a geometry file.
    IF (Scale.EQ.1) THEN
        CALL GetInputDir(InputDir,ierr)
        CALL GetInputFileTplgy(InputFileTplgy,ierr)
    END IF

    ! It obtains a temporal vector to store topology identificatiors in application ordering
    CALL DMCreateGlobalVector(DataMngr,TmpTplgy,ierr)

    ! It quantifies Dirichlet, and Cauchy boundary condition on global processors.
    BCLenL(:)=0
    BCLenG(:)=0

    ! It gets the global size from the geometry data manager.
    CALL DMDAGetInfo(DataMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
        & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,ierr)

    ! It tag every processor from 0 to n-1
    CALL MPI_Comm_rank(Comm,process,ierr)

    ! It obtains the global size of the domain on the first processor.
    IF (process.EQ.0) THEN
        Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFileTplgy))
        OPEN(u,FILE=TRIM(Route),STATUS='OLD',ACTION='READ')
        DO i=1,widthG(1)*widthG(2)*widthG(3)
            ValI=0
            READ(u,*)ValI
            IF (ValI.EQ.2) BCLenL(1)=BCLenL(1)+1 ! Dirichlet
            IF (ValI.EQ.3) BCLenL(2)=BCLenL(2)+1 ! Source
            IF (ValI.EQ.4) BCLenL(3)=BCLenL(3)+1 ! Cauchy
            ValR=REAL(ValI)
            CALL VecSetValue(TmpTplgy,i-1,ValR,INSERT_VALUES,ierr)
        END DO
        CLOSE(u)
    END IF

    CALL VecAssemblyBegin(TmpTplgy,ierr)
    CALL VecAssemblyEnd(TmpTplgy,ierr)

    ! Changing temporal topology identificators from application ordering to PETSc ordering
    CALL VecApplicationToPetsc(DataMngr,TmpTplgy,ierr)
    ViewName="ANISOFLOW_Tplgy"
    EventName="GetTopology"
    CALL ViewTopology(TmpTplgy,ViewName,EventName,ierr)

    ! Saving topology identificators as local vector in PETSc ordering
    CALL DMCreateLocalVector(DataMngr,Tplgy,ierr)

    CALL DMGlobalToLocalBegin(DataMngr,TmpTplgy,INSERT_VALUES,  &
        & Tplgy,ierr)
    CALL DMGlobalToLocalEnd(DataMngr,TmpTplgy,INSERT_VALUES,    &
        & Tplgy,ierr)

    CALL VecDestroy(TmpTplgy,ierr)

    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,"[GetGeometry Event] Topology Identifiers were satisfactorily created\n",ierr)

    CALL MPI_ALLREDUCE(BCLenL,BCLenG,3,MPI_INT,MPI_SUM,Comm,ierr)

    CALL GetTopologyBC(Comm,DataMngr,Tplgy,BCLenG,DirichIS,SourceIS,CauchyIS,ierr)

END SUBROUTINE GetTopology_1

! border of first layer dirichlet

SUBROUTINE GetTopology_2(Comm,DataMngr,Tplgy,DirichIS,SourceIS,CauchyIS,Scale,ierr)

    USE ANISOFLOW_Interface, ONLY : GetVerbose
    USE ANISOFLOW_View, ONLY : ViewTopology

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscdmda.h90>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    MPI_Comm,INTENT(IN)             :: Comm
    PetscInt,INTENT(IN)             :: Scale ! Here it do nothing, but it's keeped just to maintain syntax
    DM,INTENT(IN)                   :: DataMngr
    Vec,INTENT(OUT)                 :: Tplgy
    IS,INTENT(OUT)                  :: DirichIS,SourceIS,CauchyIS

    PetscReal,POINTER               :: TmpTplgyArray(:,:,:)
    PetscReal                       :: ValR
    PetscInt                        :: i,j,k,widthL(3),widthG(3),corn(3),BCLenL(3),BCLenG(3)
    Vec                             :: TmpTplgy
    PetscBool                       :: Verbose
    CHARACTER(LEN=200)              :: ViewName,EventName

    CALL GetVerbose(Verbose,ierr)

    ! It obtains a temporal Fortran array where will be filled each topology identifier.
    CALL DMCreateGlobalVector(DataMngr,TmpTplgy,ierr)
    CALL DMDAVecGetArrayF90(DataMngr,TmpTplgy,TmpTplgyArray,ierr)

    ! It quantifies Dirichlet, Neumman, and Cauchy boundary condition on local and global processor.
    BCLenL(:)=0
    BCLenG(:)=0

    ! It fills the temporal Fortran array.

    ! It gets the global size from the geometry data manager.
    CALL DMDAGetInfo(DataMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
        & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,ierr)

    CALL DMDAGetCorners(DataMngr,corn(1),corn(2),corn(3),widthL(1),  &
        & widthL(2),widthL(3),ierr)

    ! It assigns active identifier to every cell.
    ValR=1.0
    TmpTplgyArray(:,:,:)=ValR

    ! It changes the identifier depending on the position.
    DO k=corn(3),corn(3)+widthL(3)-1
        DO j=corn(2),corn(2)+widthL(2)-1
            DO i=corn(1),corn(1)+widthL(1)-1
                ValR=1.0
                ! Dirichlet on the border of the top layer.
                ! Active in nother case.
                IF ((k.EQ.0).AND.((i.EQ.0).OR.(i.EQ.(widthG(1)-1)).OR.     &
                & (j.EQ.0).OR.(j.EQ.(widthG(2)-1)))) THEN
                    ! Dirichlet
                    ValR=2.0
                    BCLenL(1)=BCLenL(1)+1
                END IF
                TmpTplgyArray(i,j,k)=ValR
            END DO
        END DO
    END DO
    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,"[GetGeometry Event] WARNING: Topology File wasn't provided. Default topology used is ...\n",ierr)

    ViewName="ANISOFLOW_Tplgy"
    EventName="GetTopology"
    CALL ViewTopology(TmpTplgy,ViewName,EventName,ierr)

    ! It moves the temporal Fortran array to a petsc vector, Tplgy, stored in Gmtry.
    CALL DMDAVecRestoreArrayF90(DataMngr,TmpTplgy,TmpTplgyArray,ierr)
    CALL DMCreateLocalVector(DataMngr,Tplgy,ierr)

    CALL DMGlobalToLocalBegin(DataMngr,TmpTplgy,INSERT_VALUES,  &
        & Tplgy,ierr)
    CALL DMGlobalToLocalEnd(DataMngr,TmpTplgy,INSERT_VALUES,    &
        & Tplgy,ierr)

    CALL VecDestroy(TmpTplgy,ierr)

    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,"[GetGeometry Event] Topology Identifiers were satisfactorily created\n",ierr)

    CALL MPI_ALLREDUCE(BCLenL,BCLenG,3,MPI_INT,MPI_SUM,Comm,ierr)

    CALL GetTopologyBC(Comm,DataMngr,Tplgy,BCLenG,DirichIS,SourceIS,CauchyIS,ierr)

END SUBROUTINE GetTopology_2

! bourders dirichlet

SUBROUTINE GetTopology_3(Comm,DataMngr,Tplgy,DirichIS,SourceIS,CauchyIS,Scale,ierr)

    USE ANISOFLOW_Interface, ONLY : GetVerbose
    USE ANISOFLOW_View, ONLY : ViewTopology

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscdmda.h90>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    MPI_Comm,INTENT(IN)             :: Comm
    PetscInt,INTENT(IN)             :: Scale ! Here it do nothing, but it's keeped just to maintain syntax
    DM,INTENT(IN)                   :: DataMngr
    Vec,INTENT(OUT)                 :: Tplgy
    IS,INTENT(OUT)                  :: DirichIS,SourceIS,CauchyIS

    PetscReal,POINTER               :: TmpTplgyArray(:,:,:)
    PetscReal                       :: ValR
    PetscInt                        :: i,j,k,widthL(3),widthG(3),corn(3),BCLenL(3),BCLenG(3)
    Vec                             :: TmpTplgy
    PetscBool                       :: Verbose
    CHARACTER(LEN=200)              :: ViewName,EventName

    CALL GetVerbose(Verbose,ierr)

    ! It obtains a temporal Fortran array where will be filled each topology identifier.
    CALL DMCreateGlobalVector(DataMngr,TmpTplgy,ierr)
    CALL DMDAVecGetArrayF90(DataMngr,TmpTplgy,TmpTplgyArray,ierr)

    ! It quantifies Dirichlet, Neumman, and Cauchy boundary condition on local and global processor.
    BCLenL(:)=0
    BCLenG(:)=0

    ! It fills the temporal Fortran array.

    ! It gets the global size from the geometry data manager.
    CALL DMDAGetInfo(DataMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
        & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,ierr)

    CALL DMDAGetCorners(DataMngr,corn(1),corn(2),corn(3),widthL(1),  &
        & widthL(2),widthL(3),ierr)

    ! It assigns active identifier to every cell.
    ValR=1.0
    TmpTplgyArray(:,:,:)=ValR

    ! It changes the identifier depending on the position.
    DO k=corn(3),corn(3)+widthL(3)-1
        DO j=corn(2),corn(2)+widthL(2)-1
            DO i=corn(1),corn(1)+widthL(1)-1
                ValR=1.0
                ! Dirichlet on the border
                ! Active in nother case.
                IF ((k.EQ.0).OR.(k.EQ.(widthG(3)-1)).OR. &
                  & (j.EQ.0).OR.(j.EQ.(widthG(2)-1)).OR. &
                  & (i.EQ.0).OR.(i.EQ.(widthG(1)-1))) THEN
                    ! Dirichlet
                    ValR=2.0
                    BCLenL(1)=BCLenL(1)+1
                END IF
                TmpTplgyArray(i,j,k)=ValR
            END DO
        END DO
    END DO
    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,"[GetGeometry Event] WARNING: Topology File wasn't provided. Default topology used is ...\n",ierr)

    ViewName="ANISOFLOW_Tplgy"
    EventName="GetTopology"
    CALL ViewTopology(TmpTplgy,ViewName,EventName,ierr)

    ! It moves the temporal Fortran array to a petsc vector, Tplgy, stored in Gmtry.
    CALL DMDAVecRestoreArrayF90(DataMngr,TmpTplgy,TmpTplgyArray,ierr)
    CALL DMCreateLocalVector(DataMngr,Tplgy,ierr)

    CALL DMGlobalToLocalBegin(DataMngr,TmpTplgy,INSERT_VALUES,  &
        & Tplgy,ierr)
    CALL DMGlobalToLocalEnd(DataMngr,TmpTplgy,INSERT_VALUES,    &
        & Tplgy,ierr)

    CALL VecDestroy(TmpTplgy,ierr)

    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,"[GetGeometry Event] Topology Identifiers were satisfactorily created\n",ierr)

    CALL MPI_ALLREDUCE(BCLenL,BCLenG,3,MPI_INT,MPI_SUM,Comm,ierr)

    CALL GetTopologyBC(Comm,DataMngr,Tplgy,BCLenG,DirichIS,SourceIS,CauchyIS,ierr)

END SUBROUTINE GetTopology_3

 !  - GetTopologyBC: It's a routine that creates the Boundary Condition Index Sets from Tplgy vector.
 !    > IN: DataMngr, Tplgy, BCLenG.
 !      + DataMngr: It's a DMDA PETSc structure that has stored the information related with a 
 !                 regular rectangular grid.
 !      + Tplgy: It's a PETSc vector type produced by DataMngr that contains a topology identifier
 !               in each cell.
 !          0: Inactive cell.
 !          1: Active cell.
 !          2: Dirichlet boundary condition cell.
 !          3: Source Q.
 !          4: Cauchy boundary condition
 !      + BCLenG: It's an Array that contains the global length of each Boundary Condition in the 
 !                following order: DirichIS, NeummanIS, CauchyIS, SourceIS.
 !    > OUT: DirichIS, NeummanIS, CauchyIS, SourceIS, ierr.
 !      + DirichIS: It's a PETSc index set that has a map between Dirichlet information and vecs 
 !                  produced by DataMngr.
 !      + SourceIS: It's a PETSc index set that has a map between Source information and vecs 
 !                   roduced by DataMngr.
 !      + CauchyIS: It's a PETSc index set that has a map between Cauchy information and vecs 
 !                  produced by DataMngr.
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.

SUBROUTINE GetTopologyBC(Comm,DataMngr,Tplgy,BCLenG,DirichIS,SourceIS,CauchyIS,ierr)

    USE ANISOFLOW_Interface, ONLY : GetVerbose

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscis.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    MPI_Comm,INTENT(IN)             :: Comm
    DM,INTENT(IN)                   :: DataMngr
    Vec,INTENT(IN)                  :: Tplgy
    PetscInt,INTENT(IN)             :: BCLenG(3)
    IS,INTENT(OUT)                  :: DirichIS,SourceIS,CauchyIS

    PetscInt,ALLOCATABLE            :: IndexDirich(:),IndexSource(:),IndexCauchy(:)
    PetscInt                        :: i,Low,High,Size,ValI,Count(3),SumaL(3),SumaG(3)
    PetscReal,POINTER               :: TmpTplgyArray(:)
    Vec                             :: TmpTplgy
    PetscBool                       :: Verbose

    CALL GetVerbose(Verbose,ierr)

    ALLOCATE(IndexDirich(BCLenG(1)))
    ALLOCATE(IndexSource(BCLenG(2)))
    ALLOCATE(IndexCauchy(BCLenG(3)))

    IndexDirich(:)=0
    IndexSource(:)=0
    IndexCauchy(:)=0

    CALL DMCreateGlobalVector(DataMngr,TmpTplgy,ierr)

    CALL DMLocalToGlobalBegin(DataMngr,Tplgy,INSERT_VALUES,  &
        & TmpTplgy,ierr)
    CALL DMLocalToGlobalEnd(DataMngr,Tplgy,INSERT_VALUES,    &
        & TmpTplgy,ierr)

    CALL VecGetOwnershipRange(TmpTplgy,Low,High,ierr)
    CALL VecGetSize(TmpTplgy,Size,ierr)

    CALL VecGetArrayReadF90(TmpTplgy,TmpTplgyArray,ierr)

    Count(:)=1

    DO i=0,Size
        SumaL(:)=0
        SumaG(:)=0
        IF ((i.GE.Low).AND.(i.LT.High)) THEN
            ValI=INT(TmpTplgyArray(i-Low+1))
            IF (ValI.EQ.2) THEN
                IndexDirich(Count(1))=i
                SumaL(1)=1
            ELSEIF (ValI.EQ.3) THEN
                IndexSource(Count(2))=i
                SumaL(2)=1
            ELSEIF (ValI.EQ.4) THEN
                IndexCauchy(Count(3))=i
                SumaL(3)=1
            END IF
        END IF
        CALL MPI_ALLREDUCE(SumaL,SumaG,3,MPI_INT,MPI_SUM,Comm,ierr)
        Count(:)=Count(:)+SumaG(:)
    END DO

    CALL VecRestoreArrayReadF90(TmpTplgy,TmpTplgyArray,ierr)

    CALL DMGlobalToLocalBegin(DataMngr,TmpTplgy,INSERT_VALUES,  &
        & Tplgy,ierr)
    CALL DMGlobalToLocalEnd(DataMngr,TmpTplgy,INSERT_VALUES,    &
        & Tplgy,ierr)

    CALL VecDestroy(TmpTplgy,ierr)

    CALL ISCreateGeneral(Comm,BCLenG(1),IndexDirich,                 &
        & PETSC_COPY_VALUES,DirichIS,ierr)
    CALL ISCreateGeneral(Comm,BCLenG(2),IndexSource,               &
        & PETSC_COPY_VALUES,SourceIS,ierr)
    CALL ISCreateGeneral(Comm,BCLenG(3),IndexCauchy,               &
        & PETSC_COPY_VALUES,CauchyIS,ierr)


    DEALLOCATE(IndexDirich)
    DEALLOCATE(IndexSource)
    DEALLOCATE(IndexCauchy)

    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,"[GetGeometry Event] Boundary condition index sets of topology identifiers were satisfactorily created\n",ierr)

END SUBROUTINE GetTopologyBC

SUBROUTINE GetLocalTopology(Gmtry,Ppt,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,Property,RunOptionsVar
    USE ANISOFLOW_Interface, ONLY : GetRunOptions,GetVerbose

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscdmda.h90>

    TYPE(Geometry),INTENT(IN)           :: Gmtry
    TYPE(Property),INTENT(INOUT)        :: Ppt
    PetscErrorCode,INTENT(INOUT)        :: ierr

    PetscReal,POINTER                   :: TmpTplgyArray(:,:,:),xArray(:),yArray(:),zArray(:)
    PetscInt                            :: i,j,k,xSize,ySize,zSize
    TYPE(RunOptionsVar)                 :: RunOptions

    i=Ppt%Pstn%i
    j=Ppt%Pstn%j
    k=Ppt%Pstn%k

    CALL VecGetSize(Gmtry%x,xSize,ierr)
    CALL VecGetArrayReadF90(Gmtry%x,xArray,ierr)
    Ppt%dxB=xArray(i+1)-xArray(i)
    Ppt%dx =xArray(i+2)-xArray(i+1)
    Ppt%dxF=xArray(i+3)-xArray(i+2)
    CALL VecRestoreArrayReadF90(Gmtry%x,xArray,ierr)

    CALL VecGetSize(Gmtry%y,ySize,ierr)
    CALL VecGetArrayReadF90(Gmtry%y,yArray,ierr)
    Ppt%dyB=yArray(j+1)-yArray(j)
    Ppt%dy =yArray(j+2)-yArray(j+1)
    Ppt%dyF=yArray(j+3)-yArray(j+2)
    CALL VecRestoreArrayReadF90(Gmtry%y,yArray,ierr)

    CALL VecGetSize(Gmtry%z,zSize,ierr)
    CALL VecGetArrayReadF90(Gmtry%z,zArray,ierr)
    Ppt%dyB=zArray(k+1)-zArray(k)
    Ppt%dy =zArray(k+2)-zArray(k+1)
    Ppt%dyF=zArray(k+3)-zArray(k+2)
    CALL VecRestoreArrayReadF90(Gmtry%z,zArray,ierr)

    ! Revisar los diferenciales de x,y,z que se usan como dividendo, si los valores son cercanos a cero puede haber problemas!!!!

    CALL GetRunOptions(RunOptions,ierr)

    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,Gmtry%Tplgy,TmpTplgyArray,ierr)

    IF (RunOptions%Scheme.EQ.1) THEN
        ALLOCATE(Ppt%StnclTplgy(7))
        Ppt%StnclTplgy(:)=0
        Ppt%StnclTplgy(1)=INT(TmpTplgyArray(i  ,j  ,k-1))
        Ppt%StnclTplgy(2)=INT(TmpTplgyArray(i  ,j-1,k  ))
        Ppt%StnclTplgy(3)=INT(TmpTplgyArray(i-1,j  ,k  ))
        Ppt%StnclTplgy(4)=INT(TmpTplgyArray(i  ,j  ,k  ))
        Ppt%StnclTplgy(5)=INT(TmpTplgyArray(i+1,j  ,k  ))
        Ppt%StnclTplgy(6)=INT(TmpTplgyArray(i  ,j+1,k  ))
        Ppt%StnclTplgy(7)=INT(TmpTplgyArray(i  ,j  ,k+1))
    ELSEIF (RunOptions%Scheme.EQ.2) THEN
        ALLOCATE(Ppt%StnclTplgy(19))
        Ppt%StnclTplgy(:)=0
        Ppt%StnclTplgy(1)= INT(TmpTplgyArray(i  ,j-1,k-1))
        Ppt%StnclTplgy(2)= INT(TmpTplgyArray(i-1,j  ,k-1))
        Ppt%StnclTplgy(3)= INT(TmpTplgyArray(i  ,j  ,k-1))
        Ppt%StnclTplgy(4)= INT(TmpTplgyArray(i+1,j  ,k-1))
        Ppt%StnclTplgy(5)= INT(TmpTplgyArray(i  ,j+1,k-1))
        Ppt%StnclTplgy(6)= INT(TmpTplgyArray(i-1,j-1,k  ))
        Ppt%StnclTplgy(7)= INT(TmpTplgyArray(i  ,j-1,k  ))
        Ppt%StnclTplgy(8)= INT(TmpTplgyArray(i+1,j-1,k  ))
        Ppt%StnclTplgy(9)= INT(TmpTplgyArray(i-1,j  ,k  ))
        Ppt%StnclTplgy(10)=INT(TmpTplgyArray(i  ,j  ,k  ))
        Ppt%StnclTplgy(11)=INT(TmpTplgyArray(i+1,j  ,k  ))
        Ppt%StnclTplgy(12)=INT(TmpTplgyArray(i-1,j+1,k  ))
        Ppt%StnclTplgy(13)=INT(TmpTplgyArray(i  ,j+1,k  ))
        Ppt%StnclTplgy(14)=INT(TmpTplgyArray(i+1,j+1,k  ))
        Ppt%StnclTplgy(15)=INT(TmpTplgyArray(i  ,j-1,k+1))
        Ppt%StnclTplgy(16)=INT(TmpTplgyArray(i-1,j  ,k+1))
        Ppt%StnclTplgy(17)=INT(TmpTplgyArray(i  ,j  ,k+1))
        Ppt%StnclTplgy(18)=INT(TmpTplgyArray(i+1,j  ,k+1))
        Ppt%StnclTplgy(19)=INT(TmpTplgyArray(i  ,j+1,k+1))
    ELSE
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "[ERROR] Run_options_scheme command must be an integer between 1 and 2\n",ierr)
        STOP
    END IF
    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,Gmtry%Tplgy,TmpTplgyArray,ierr)

END SUBROUTINE GetLocalTopology

SUBROUTINE DestroyGeometry(Gmtry,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry
    USE ANISOFLOW_Interface, ONLY : GetVerbose

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscis.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    TYPE(Geometry),INTENT(INOUT)    :: Gmtry

    PetscBool                       :: Verbose
    CHARACTER(LEN=200)              :: EventName,ClassName
    PetscLogEvent                   :: Event
    PetscClassId                    :: ClassID
    PetscLogDouble                  :: EventFlops

    ClassName="Geometry"
    CALL PetscClassIdRegister(ClassName,ClassID,ierr)
    EventName="DestroyGeometry"
    CALL PetscLogEventRegister(EventName,ClassID,Event,ierr)
    CALL PetscLogEventBegin(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Inizialited\n",ierr)

    CALL DMDestroy(Gmtry%DataMngr,ierr)
    CALL VecDestroy(Gmtry%Tplgy,ierr)
    CALL VecDestroy(Gmtry%x,ierr)
    CALL VecDestroy(Gmtry%y,ierr)
    CALL VecDestroy(Gmtry%z,ierr)
    CALL ISDestroy(Gmtry%DirichIS,ierr)
    CALL ISDestroy(Gmtry%SourceIS,ierr)
    CALL ISDestroy(Gmtry%CauchyIS,ierr)

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Finalized\n",ierr)
    
    CALL PetscLogFlops(EventFlops,ierr)
    CALL PetscLogEventEnd(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)


END SUBROUTINE DestroyGeometry

SUBROUTINE VecApplicationToPetsc(DataMngr,AppVec,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscis.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscao.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    DM,INTENT(IN)                   :: DataMngr
    Vec,INTENT(INOUT)               :: AppVec ! it outs as PetscVec

    AO                              :: AppOrd
    PetscInt                        :: RangeLow,RangeHigh,i
    IS                              :: PetscIS,AppIS
    VecScatter                      :: Scatter
    Vec                             :: PetscVec

    ! It obtains a temporal vector to store in PETSc ordering
    CALL DMCreateGlobalVector(DataMngr,PetscVec,ierr)
    ! Changing temporal vector from application ordering to PETSc ordering
    CALL DMDAGetAO(DataMngr,AppOrd,ierr)
    CALL VecGetOwnershipRange(PetscVec,RangeLow,RangeHigh,ierr)
    CALL ISCreateGeneral(PETSC_COMM_WORLD,RangeHigh-RangeLow,(/(i,i=RangeLow,RangeHigh)/),PETSC_COPY_VALUES,PetscIS,ierr)
    CALL ISDuplicate(PetscIS,AppIS,ierr)
    CALL ISCopy(PetscIS,AppIS,ierr)
    CALL AOPetscToApplicationIS(AppOrd,AppIS,ierr)
    CALL VecScatterCreate(AppVec,AppIS,PetscVec,PetscIS,Scatter,ierr)
    CALL VecScatterBegin(Scatter,AppVec,PetscVec,INSERT_VALUES,SCATTER_FORWARD,ierr)
    CALL VecScatterEnd(Scatter,AppVec,PetscVec,INSERT_VALUES,SCATTER_FORWARD,ierr)

    CALL VecScatterDestroy(Scatter,ierr)
    CALL ISDestroy(AppIS,ierr)
    CALL ISDestroy(PetscIS,ierr)
    CALL VecCopy(PetscVec,AppVec,ierr)
    CALL VecDestroy(PetscVec,ierr)

END SUBROUTINE VecApplicationToPetsc

END MODULE ANISOFLOW_Geometry