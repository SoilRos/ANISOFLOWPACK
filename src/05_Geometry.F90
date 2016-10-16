MODULE ANISOFLOW_Geometry

! ANISOFLOW_Geometry it's a module that contains routines to manage 
! geometry variables.

    IMPLICIT NONE

CONTAINS

 !  - GetGeometry: It's a routine that fills a Geometry data structure 
 !                 with input files provided by
 !                 the user.
 !    > OUT: Gmtry, ierr.
 !      + Gmtry: It's a Geometry data structure filled with input 
 !               files provided by the user.
 !      + ierr: It's an integer that indicates whether an error has 
 !              occurred during the call.
 !    > NOTES: The Geometry is filled in three stages: to create a 
 !             Data Manager (DataMngr) to control information related 
 !             to a regular rectangular grid, a grid indformation and, 
 !             a Topology that describes the geometry with 
 !             identifiers.

SUBROUTINE GetGeometry(Comm,Gmtry,ierr)

    USE ANISOFLOW_Interface,    ONLY : GetVerbose
    USE ANISOFLOW_Types,        ONLY : Geometry

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    MPI_Comm,INTENT(IN)             :: Comm
    TYPE(Geometry),INTENT(OUT)      :: Gmtry

    PetscBool                       :: Verbose
    CHARACTER(LEN=200)              :: EventName,ClassName
    PetscLogEvent                   :: Event
    PetscClassId                    :: ClassID
    PetscLogDouble                  :: EventFlops=0.d0

    ClassName="Geometry"
    CALL PetscClassIdRegister(ClassName,ClassID,ierr)
    EventName="GetGeometry"
    CALL PetscLogEventRegister(EventName,ClassID,Event,ierr)
    CALL PetscLogEventBegin(Event,PETSC_NULL_OBJECT,                 &
        & PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,"["//             &
        & ADJUSTL(TRIM(EventName))//" Event] Inizialited\n",ierr)
    
    CALL GetDataMngr(Comm,Gmtry%Scale,Gmtry%DataMngr,ierr)
    CALL GetGrid(Comm,Gmtry%DataMngr,Gmtry%Scale,Gmtry%x,Gmtry%y,    &
        & Gmtry%z,ierr)
    CALL GetTopology(Comm,Gmtry%DataMngr,Gmtry%Scale,Gmtry%Tplgy,    &
        & Gmtry%SizeTplgy,Gmtry%InactiveIS,ierr)

    ! Storing a permanet topology on pTplgy
    CALL VecDuplicate(Gmtry%Tplgy,Gmtry%pTplgy,ierr)
    CALL VecCopy(Gmtry%Tplgy,Gmtry%pTplgy,ierr)

    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,"["//             &
        & ADJUSTL(TRIM(EventName))//" Event] Finalized\n",ierr)
    
    CALL PetscGetFlops(EventFlops,ierr)
    CALL PetscLogFlops(EventFlops,ierr)
    CALL PetscLogEventEnd(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
        & PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

END SUBROUTINE GetGeometry

 !  - GetDataMngr: It's a routine that creates and fills the 
 !                 information related with a regular rectangular 
 !                 grid.
 !    > OUT: DataMngr, ierr.
 !      + DataMngr: It's a DMDA PETSc structure that stores the 
 !                  information related with a regular rectangular 
 !                  grid.
 !      + ierr: It's an integer that indicates whether an error has 
 !              occurred during the call.
 !    > NOTES: The DataMngr takes the size of the domain and assign a 
 !             subdomain to each processor.

SUBROUTINE GetDataMngr(Comm,Scale,DataMngr,ierr)

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

    ! InputType defines the type of file that is provided.
    !   1 and 2: Defined by Perez.

    IF ((InputType%Gmtry.EQ.1).OR.(InputType%Gmtry.EQ.2)) THEN
        CALL GetDataMngr_1(Comm,Scale,DataMngr,ierr)
    ELSE
        CALL PetscSynchronizedPrintf(Comm, &
            & "[ERROR] Geometry InputType wrong\n",ierr)
        STOP
    END IF

END SUBROUTINE GetDataMngr

SUBROUTINE GetDataMngr_1(Comm,Scale,DataMngr,ierr)

    USE ANISOFLOW_Types,        ONLY : RunOptionsVar
    USE ANISOFLOW_Interface,    ONLY : GetInputDir,GetInputFileGmtry,&
                                     & GetRunOptions,GetVerbose

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
    CHARACTER(LEN=200)              :: InputDir,InputFileGmtry,Route,&
                                     & aux
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
        READ(u,*)aux,widthG(1)
        READ(u,*)aux,widthG(2)
        READ(u,*)aux,widthG(3)
        CLOSE(u)
    END IF

    ! It broadcasts the global size to other processors.
    CALL MPI_Bcast(widthG,3,MPI_INT, 0, Comm,ierr)

    ! It decides the stencil shape depending on the scheme used.
    IF (RunOptions%Scheme.EQ.1) THEN
        Stencil=DMDA_STENCIL_STAR
    ELSEIF ((RunOptions%Scheme.EQ.2).OR.(RunOptions%Scheme.EQ.3)) THEN
        Stencil=DMDA_STENCIL_BOX
    ELSE
        CALL PetscSynchronizedPrintf(Comm,                           &
            & "[ERROR] Run_options_scheme command must be an integ"//&
            & "er between 1 and 3.\n",ierr)
        STOP
    END IF

    ! It creates the Data Manager to Distributed Arrays by the 
    ! information provided.
    IF (widthG(3).NE.1) THEN
        CALL DMDACreate3d(Comm,DM_BOUNDARY_GHOSTED,                  &
            & DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,Stencil,       &
            & widthG(1),widthG(2),widthG(3),PETSC_DECIDE,            &
            & PETSC_DECIDE,PETSC_DECIDE,1,1,PETSC_NULL_INTEGER,      &
            & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,DataMngr,ierr)
    ELSE
        CALL DMDACreate2d(Comm,DM_BOUNDARY_GHOSTED,                  &
            & DM_BOUNDARY_GHOSTED,Stencil,widthG(1),widthG(2),       &
            & PETSC_DECIDE,PETSC_DECIDE,1,1,PETSC_NULL_INTEGER,      &
            & PETSC_NULL_INTEGER,DataMngr,ierr)
    END IF

    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,                  &
        & "[GetGeometry Event] Data Manager was satisfactorily cre"//&
        & "ated.\n",ierr)

END SUBROUTINE GetDataMngr_1

 !  - GetGrid: It's a routine that creates and fills the information 
 !             related with a a coordinates of a regular rectangular 
 !             grid.
 !    > IN: DataMngr
 !      + DataMngr: It's a DMDA PETSc structure that has stored the 
 !                  information related with a regular rectangular 
 !                  grid.
 !    > OUT: x, y, z, ierr.
 !      + x,y,z: It's an Array that has de coordinates of the grid on 
 !               each direction. It's stored on every processor 
 !      + ierr: It's an integer that indicates whether an error has 
 !              occurred during the call.

SUBROUTINE GetGrid(Comm,DataMngr,Scale,x,y,z,ierr)

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
    !   2: Defined by Perez.

    IF (InputType%Gmtry.EQ.1) THEN
        CALL GetGrid_1(Comm,DataMngr,Scale,x,y,z,ierr)
    ELSE IF (InputType%Gmtry.EQ.2) THEN
        CALL GetGrid_2(Comm,DataMngr,Scale,x,y,z,ierr)
    ELSE        
        CALL PetscSynchronizedPrintf(Comm, &
            & "[ERROR] Geometry InputType wrong\n",ierr)
        STOP
    END IF

END SUBROUTINE GetGrid

 !  - GetGrid_1: It's a routine that creates and fills the information 
 !               related with a a coordinates of
 !               a regular rectangular grid when InputType%Gmtry=1
 !    > IN: DataMngr
 !      + DataMngr: It's a DMDA PETSc structure that has stored the 
 !                  information related with a 
 !                 regular rectangular grid.
 !    > OUT: x, y, z, ierr.
 !      + x,y,z: It's an Array that has de coordinates of the grid on 
 !               each direction. It's stored on 
 !               every processor 
 !      + ierr: It's an integer that indicates whether an error has 
 !              occurred during the call.

SUBROUTINE GetGrid_1(Comm,DataMngr,Scale,x,y,z,ierr)

    USE ANISOFLOW_Interface, ONLY : GetInputDir,GetInputFileGmtry,   &
                                  & GetVerbose

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
        & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,           &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,  &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,  &
        & PETSC_NULL_INTEGER,ierr)

    ! It create the vectors
    size=widthG(1)+1
    CALL VecCreateSeq(PETSC_COMM_SELF,size,x,ierr)

    size=widthG(2)+1
    CALL VecCreateSeq(PETSC_COMM_SELF,size,y,ierr)

    size=widthG(3)+1
    CALL VecCreateSeq(PETSC_COMM_SELF,size,z,ierr)

    
    DO i=0,widthG(1)
        Value=REAL(i)
        CALL VecSetValue(x,i,Value,INSERT_VALUES,ierr)
    END DO
    DO i=0,widthG(2)
        Value=REAL(i)
        CALL VecSetValue(y,i,Value,INSERT_VALUES,ierr)
    END DO
    DO i=0,widthG(3)
        Value=REAL(i)
        CALL VecSetValue(z,i,Value,INSERT_VALUES,ierr)
    END DO
    ! Default DX=DY=DZ=1.0
    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,                  &
        & "[GetGeometry Event] WARNING: Grid System wasn't provide"//&
        & "d. To provide one you have to use -Input_type_gmtry 2 o"//&
        & "therwise a default grid used is DX=DY=DZ=1.0.\n",ierr)

    CALL VecAssemblyBegin(x,ierr)
    CALL VecAssemblyEnd(x,ierr)

    CALL VecAssemblyBegin(y,ierr)
    CALL VecAssemblyEnd(y,ierr)

    CALL VecAssemblyBegin(z,ierr)
    CALL VecAssemblyEnd(z,ierr)

    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,                  &
        & "[GetGeometry Event] Grid System was satisfactorily crea"//&
        & "ted.\n",ierr)

END SUBROUTINE GetGrid_1

SUBROUTINE GetGrid_2(Comm,DataMngr,Scale,x,y,z,ierr)

    USE ANISOFLOW_Interface, ONLY : GetInputDir,GetInputFileGmtry,   &
                                  & GetVerbose

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    MPI_Comm,INTENT(IN)                 :: Comm
    PetscInt,INTENT(IN)                 :: Scale
    DM,INTENT(IN)                       :: DataMngr
    Vec,INTENT(OUT)                     :: x,y,z

    PetscMPIInt                         :: process
    CHARACTER(LEN=200)                  :: InputDir,InputFileGmtry,  &
                                         & Route
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
        & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,           &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,  &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,  &
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
            READ(u,*)Value
        END IF
        ! It broadcasts the global size to other processors.
        CALL MPI_Bcast(Value,1,MPI_DOUBLE, 0, Comm,ierr)
        CALL VecSetValue(x,i,Value,INSERT_VALUES,ierr)
    END DO
    
    DO i=0,widthG(2)
        IF (process.EQ.0) THEN
            READ(u,*)Value
        END IF
        ! It broadcasts the global size to other processors.
        CALL MPI_Bcast(Value,1,MPI_DOUBLE, 0, Comm,ierr)
        CALL VecSetValue(y,i,Value,INSERT_VALUES,ierr)
    END DO

    DO i=0,widthG(3)
        IF (process.EQ.0) THEN
            READ(u,*)Value
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

    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,                  &
        & "[GetGeometry Event] Grid System was satisfactorily crea"//&
        & "ted\n",ierr)

END SUBROUTINE GetGrid_2

 !  - GetTopology: It's a routine that creates and fills the 
 !                 information related to topology. It creates a      
 !                 vector and index sets to describe the geometry. 
 !    > IN: DataMngr.
 !      + DataMngr: It's a DMDA PETSc structure that has stored the 
 !                  information related with a regular rectangular 
 !                  grid.
 !    > OUT: Tplgy, DirichIS, NeummanIS, CauchyIS, SourceIS, ierr.
 !      + Tplgy: It's a PETSc vector type produced by DataMngr that 
 !               contains a topology identifier in each cell.
 !          0: Inactive cell.
 !          1: Active cell.
 !          2: Dirichlet boundary condition cell.
 !          3: Source Q.
 !          4: Cauchy boundary condition
 !      + SizeTplgy: It's a array of integers that quantify the 
 !                   indentifiers of each type decribed above where 
 !                   inactive cells quantifier doesn't exist.
 !          1: Active cell quantifier.
 !          2: Dirichlet boundary condition cell quantifier.
 !          3: Source in cell quantifier.
 !          4: Cauchy boundary condition quantifier.
 !      + ierr: It's an integer that indicates whether an error has 
 !              occurred during the call.

SUBROUTINE GetTopology(Comm,DataMngr,Scale,Tplgy,SizeTplgy,          &
    & InactiveIS,ierr)

    USE ANISOFLOW_Types, ONLY : InputTypeVar
    USE ANISOFLOW_Interface, ONLY : GetInputType,GetVerbose

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    MPI_Comm,INTENT(IN)             :: Comm
    PetscInt,INTENT(IN)             :: Scale
    DM,INTENT(IN)                   :: DataMngr
    Vec,INTENT(OUT)                 :: Tplgy
    PetscInt,INTENT(OUT)            :: SizeTplgy(4)
    IS,INTENT(OUT)                  :: InactiveIS

    TYPE(InputTypeVar)              :: InputType

    IF (Scale.EQ.1) THEN
        CALL GetInputType(InputType,ierr)
    END IF

    ! InputType define the type of file that is provided.
    !   0: all active
    !   1: From file
    !   2: Border of first layer dirichlet
    !   3: Borders dirichlet
    IF (InputType%Tplgy.EQ.0) THEN
        CALL GetTopology_0(Comm,DataMngr,Scale,Tplgy,SizeTplgy,      &
            & InactiveIS,ierr)
    ELSEIF (InputType%Tplgy.EQ.1) THEN
        CALL GetTopology_1(Comm,DataMngr,Scale,Tplgy,SizeTplgy,      &
            & InactiveIS,ierr)
    ELSE IF (InputType%Tplgy.EQ.2) THEN
        CALL GetTopology_2(Comm,DataMngr,Scale,Tplgy,SizeTplgy,      &
            & InactiveIS,ierr)
    ELSE IF (InputType%Tplgy.EQ.3) THEN
        CALL GetTopology_3(Comm,DataMngr,Scale,Tplgy,SizeTplgy,      &
            & InactiveIS,ierr)
    ELSE
        CALL PetscSynchronizedPrintf(Comm, &
            & "[ERROR] Topology InputType wrong\n",ierr)
        STOP
    END IF

    IF ((SizeTplgy(2)+SizeTplgy(3)+SizeTplgy(4)).EQ.0) THEN
        ! Imprimir un mensaje de advertencia, no esta mal, pero obliga
        ! a que se agregen condiciones de frontera "temporales".
    END IF

END SUBROUTINE GetTopology


SUBROUTINE GetTopology_0(Comm,DataMngr,Scale,Tplgy,SizeTplgy,        & 
    & InactiveIS,ierr)

    USE ANISOFLOW_Interface, ONLY : GetVerbose,GetInputDir,          &
                                  & GetInputFileTplgy
    USE ANISOFLOW_View,      ONLY : ViewTopology

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscis.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscdmda.h90>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    MPI_Comm,INTENT(IN)             :: Comm
    PetscInt,INTENT(IN)             :: Scale
    DM,INTENT(IN)                   :: DataMngr
    Vec,INTENT(OUT)                 :: Tplgy
    PetscInt,INTENT(OUT)            :: SizeTplgy(4)
    IS,INTENT(OUT)                  :: InactiveIS

    CHARACTER(LEN=200)              :: ViewName,EventName
    PetscReal                       :: one=1.D0
    PetscInt                        :: widthG(3),SizeInactive
    PetscInt,ALLOCATABLE            :: IndexInactive(:)
    PetscBool                       :: Verbose


    CALL GetVerbose(Verbose,ierr)

    ! It quantifies Dirichlet, and Cauchy boundary condition on global
    ! processors.
    SizeInactive=0
    SizeTplgy(:)=0

    ! It gets the global size from the geometry data manager.
    CALL DMDAGetInfo(DataMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
        & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,           &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,  &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,  &
        & PETSC_NULL_INTEGER,ierr)

    ! Saving topology identificators as local vector in PETSc ordering
    CALL DMCreateLocalVector(DataMngr,Tplgy,ierr)

    CALL VecSet(Tplgy,one,ierr)

    ! Changing topology identificators from application ordering to 
    ! PETSc ordering
    ViewName="Tplgy"
    EventName="GetTopology"
    CALL ViewTopology(Tplgy,ViewName,EventName,ierr)

    SizeTplgy(1)=widthG(1)*widthG(2)*widthG(3)

    ALLOCATE(IndexInactive(SizeInactive))
    CALL MPI_Bcast(IndexInactive,SizeInactive,MPI_INT,0,Comm,ierr)
    CALL ISCreateGeneral(Comm,SizeInactive,IndexInactive,            &
        & PETSC_COPY_VALUES,InactiveIS,ierr)
    DEALLOCATE(IndexInactive)

    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,                  &
        & "[GetGeometry Event] WARNING: All domain is active.\n"     &
        & ,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,                  &
        & "[GetGeometry Event] Topology Identifiers were satisfact"//&
        & "orily created\n",ierr)

END SUBROUTINE GetTopology_0

 !  - GetTopology_1: It's a routine that creates and fills the 
 !                   information related to topology. when 
 !                   InputType%Tplgy=1. It creates a vector and index 
 !                   sets to describe the geometry. 
 !    > IN: DataMngr.
 !      + DataMngr: It's a DMDA PETSc structure that has stored the 
 !                  information related with a regular rectangular 
 !                  grid. 
 !    > OUT: Tplgy, DirichIS, NeummanIS, CauchyIS, SourceIS, ierr.
 !      + Tplgy: It's a PETSc vector type produced by DataMngr that 
 !               contains a topology identifier in each cell.
 !          0: Inactive cell.
 !          1: Active cell.
 !          2: Dirichlet boundary condition cell.
 !          3: Source Q.
 !          4: Cauchy boundary condition
 !      + SizeTplgy: It's a array of integers that quantify the 
 !                   indentifiers of each type decribed above where 
 !                   inactive cells quantifier doesn't exist. It can 
 !                   be calculated with  the rest.
 !          1: Active cell quantifier.
 !          2: Dirichlet boundary condition cell quantifier.
 !          3: Source in cell quantifier.
 !          4: Cauchy boundary condition quantifier.
 !      + ierr: It's an integer that indicates whether an error has 
 !              occurred during the call.

SUBROUTINE GetTopology_1(Comm,DataMngr,Scale,Tplgy,SizeTplgy,        & 
    & InactiveIS,ierr)

    USE ANISOFLOW_Interface, ONLY : GetVerbose,GetInputDir,          &
                                  & GetInputFileTplgy
    USE ANISOFLOW_View,      ONLY : ViewTopology

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscis.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscdmda.h90>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    MPI_Comm,INTENT(IN)             :: Comm
    PetscInt,INTENT(IN)             :: Scale
    DM,INTENT(IN)                   :: DataMngr
    Vec,INTENT(OUT)                 :: Tplgy
    PetscInt,INTENT(OUT)            :: SizeTplgy(4)
    IS,INTENT(OUT)                  :: InactiveIS

    PetscMPIInt                     :: process
    CHARACTER(LEN=200)              :: InputDir,InputFileTplgy,Route,&
                                     & ViewName,EventName
    PetscReal                       :: ValR
    PetscInt                        :: u,i,ValI,widthG(3),SizeInactive
    PetscInt,ALLOCATABLE            :: IndexInactive(:),             &
                                     & TmpIndexInactive(:)
    Vec                             :: TmpTplgy
    PetscBool                       :: Verbose

    PARAMETER(u=01)

    CALL GetVerbose(Verbose,ierr)

    ! It obtains the route to open a geometry file.
    IF (Scale.EQ.1) THEN
        CALL GetInputDir(InputDir,ierr)
        CALL GetInputFileTplgy(InputFileTplgy,ierr)
    END IF

    ! It obtains a temporal vector to store topology identificatiors 
    ! in application ordering
    CALL DMCreateGlobalVector(DataMngr,TmpTplgy,ierr)

    ! It quantifies Dirichlet, and Cauchy boundary condition on global 
    ! processors.
    SizeInactive=0
    SizeTplgy(:)=0

    ! It gets the global size from the geometry data manager.
    CALL DMDAGetInfo(DataMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
        & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,           &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,  &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,  &
        & PETSC_NULL_INTEGER,ierr)

    ! It tag every processor from 0 to n-1
    CALL MPI_Comm_rank(Comm,process,ierr)


    ! It obtains the global size of the domain on the first processor.
    IF (process.EQ.0) THEN
        Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFileTplgy))
        OPEN(u,FILE=TRIM(Route),STATUS='OLD',ACTION='READ')
        ALLOCATE(TmpIndexInactive(widthG(1)*widthG(2)*widthG(3)))
        DO i=1,widthG(1)*widthG(2)*widthG(3)
            ValI=0
            READ(u,*)ValI
            IF (ValI.EQ.0) SizeInactive=SizeInactive+1
            IF (ValI.EQ.0) TmpIndexInactive(SizeInactive)=i-1
            IF (ValI.EQ.1) SizeTplgy(1)=SizeTplgy(1)+1 ! Active
            IF (ValI.EQ.2) SizeTplgy(2)=SizeTplgy(2)+1 ! Dirichlet
            IF (ValI.EQ.3) SizeTplgy(3)=SizeTplgy(3)+1 ! Source
            IF (ValI.EQ.4) SizeTplgy(4)=SizeTplgy(4)+1 ! Cauchy
            ValR=REAL(ValI)
            CALL VecSetValue(TmpTplgy,i-1,ValR,INSERT_VALUES,ierr)
        END DO
        CLOSE(u)
    END IF

    CALL VecAssemblyBegin(TmpTplgy,ierr)
    CALL VecAssemblyEnd(TmpTplgy,ierr)

    ! Changing temporal topology identificators from application 
    ! ordering to PETSc ordering
    CALL VecApplicationToPetsc(DataMngr,TmpTplgy,ierr)
    ViewName="Tplgy"
    EventName="GetTopology"
    CALL ViewTopology(TmpTplgy,ViewName,EventName,ierr)

    ! Saving topology identificators as local vector in PETSc ordering
    CALL DMCreateLocalVector(DataMngr,Tplgy,ierr)

    CALL DMGlobalToLocalBegin(DataMngr,TmpTplgy,INSERT_VALUES,Tplgy, &
        & ierr)
    CALL DMGlobalToLocalEnd(DataMngr,TmpTplgy,INSERT_VALUES,Tplgy,ierr)

    CALL VecDestroy(TmpTplgy,ierr)

    CALL MPI_Bcast(SizeTplgy,4,MPI_INT,0,Comm,ierr)
    CALL MPI_Bcast(SizeInactive,1,MPI_INT,0,Comm,ierr)

    IF (SizeInactive.NE.(widthG(1)*widthG(2)*widthG(3)-(SizeTplgy(1)+&
        & SizeTplgy(2)+SizeTplgy(3)+SizeTplgy(4)))) THEN
        CALL PetscSynchronizedPrintf(Comm,                           &
            & "[GetGeometry Event] ERROR: Topology file is invalid"//&
            & ". Maybe due integers in the file are not between 0 "//&
            & "and 4.\n",ierr)
    END IF
    ALLOCATE(IndexInactive(SizeInactive))
    IF (process.EQ.0) IndexInactive(:)=TmpIndexInactive(1:SizeInactive)
    IF (process.EQ.0) DEALLOCATE(TmpIndexInactive)

    CALL MPI_Bcast(IndexInactive,SizeInactive,MPI_INT,0,Comm,ierr)
    CALL ISCreateGeneral(Comm,SizeInactive,IndexInactive,            &
        & PETSC_COPY_VALUES,InactiveIS,ierr)
    DEALLOCATE(IndexInactive)

    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,                  &
        & "[GetGeometry Event] Topology Identifiers were satisfact"//&
        & "orily created.\n",ierr)

END SUBROUTINE GetTopology_1

SUBROUTINE GetTopology_2(Comm,DataMngr,Scale,Tplgy,SizeTplgy,        &
    & InactiveIS,ierr)

    USE ANISOFLOW_Interface, ONLY : GetVerbose
    USE ANISOFLOW_View,      ONLY : ViewTopology

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscis.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscdmda.h90>

    ! Scale do nothing here, it's keeped just to maintain the syntax.
    PetscErrorCode,INTENT(INOUT)    :: ierr
    MPI_Comm,INTENT(IN)             :: Comm
    PetscInt,INTENT(IN)             :: Scale 
    DM,INTENT(IN)                   :: DataMngr
    Vec,INTENT(OUT)                 :: Tplgy
    PetscInt,INTENT(OUT)            :: SizeTplgy(4)
    IS,INTENT(OUT)                  :: InactiveIS

    PetscReal,POINTER               :: TmpTplgyArray(:,:,:)
    PetscReal                       :: ValR
    PetscInt                        :: i,j,k,widthL(3),widthG(3),    &
                                     & corn(3),SizeInactive
    PetscInt,ALLOCATABLE            :: IndexInactive(:)
    Vec                             :: TmpTplgy
    PetscBool                       :: Verbose
    CHARACTER(LEN=200)              :: ViewName,EventName

    CALL GetVerbose(Verbose,ierr)

    ! It obtains a temporal Fortran array where will be filled each 
    ! topology identifier.
    CALL DMCreateGlobalVector(DataMngr,TmpTplgy,ierr)
    CALL DMDAVecGetArrayF90(DataMngr,TmpTplgy,TmpTplgyArray,ierr)

    ! It quantifies Dirichlet, Neumman, and Cauchy boundary condition 
    ! on local and global processor.
    SizeTplgy(:)=0

    ! It fills the temporal Fortran array.

    ! It gets the global size from the geometry data manager.
    CALL DMDAGetInfo(DataMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
        & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,           &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,  &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,  &
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
                IF ((k.EQ.0).AND.((i.EQ.0).OR.(i.EQ.(widthG(1)-1))   &
                    &.OR.(j.EQ.0).OR.(j.EQ.(widthG(2)-1)))) THEN
                    ! Dirichlet
                    ValR=2.0
                    SizeTplgy(2)=SizeTplgy(2)+1
                END IF
                TmpTplgyArray(i,j,k)=ValR
            END DO
        END DO
    END DO
    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,                  &
        & "[GetGeometry Event] WARNING: Topology File wasn't provi"//&
        & "ded. Default topology used is ...\n",ierr)

    ViewName="Tplgy"
    EventName="GetTopology"
    CALL ViewTopology(TmpTplgy,ViewName,EventName,ierr)

    ! It moves the temporal Fortran array to a petsc vector, Tplgy, 
    ! stored in Gmtry.
    CALL DMDAVecRestoreArrayF90(DataMngr,TmpTplgy,TmpTplgyArray,ierr)
    CALL DMCreateLocalVector(DataMngr,Tplgy,ierr)

    CALL DMGlobalToLocalBegin(DataMngr,TmpTplgy,INSERT_VALUES,Tplgy, &
        & ierr)
    CALL DMGlobalToLocalEnd(DataMngr,TmpTplgy,INSERT_VALUES,Tplgy,   &
        & ierr)

    CALL VecDestroy(TmpTplgy,ierr)

    CALL MPI_Bcast(SizeTplgy,4,MPI_INT,0,Comm,ierr)
    CALL MPI_Bcast(SizeInactive,1,MPI_INT,0,Comm,ierr)
    SizeTplgy(1)=widthG(1)*widthG(2)*widthG(3)

    ALLOCATE(IndexInactive(SizeInactive))
    CALL MPI_Bcast(IndexInactive,SizeInactive,MPI_INT,0,Comm,ierr)
    CALL ISCreateGeneral(Comm,SizeInactive,IndexInactive,            &
        & PETSC_COPY_VALUES,InactiveIS,ierr)
    DEALLOCATE(IndexInactive)

    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,                  &
        & "[GetGeometry Event] Topology Identifiers were satisfact"//&
        & "orily created\n",ierr)

END SUBROUTINE GetTopology_2

SUBROUTINE GetTopology_3(Comm,DataMngr,Scale,Tplgy,SizeTplgy,        &
    & InactiveIS,ierr)

    USE ANISOFLOW_Interface, ONLY : GetVerbose
    USE ANISOFLOW_View,      ONLY : ViewTopology

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscdmda.h90>

    ! Scale do nothing here, it's keeped just to maintain the syntax.
    PetscErrorCode,INTENT(INOUT)    :: ierr
    MPI_Comm,INTENT(IN)             :: Comm
    PetscInt,INTENT(IN)             :: Scale
    DM,INTENT(IN)                   :: DataMngr
    Vec,INTENT(OUT)                 :: Tplgy
    PetscInt,INTENT(OUT)            :: SizeTplgy(4)
    IS,INTENT(OUT)                  :: InactiveIS

    PetscReal,POINTER               :: TmpTplgyArray(:,:,:)
    PetscReal                       :: ValR
    PetscInt                        :: i,j,k,widthL(3),widthG(3),    &
                                     & corn(3),SizeInactive
    PetscInt,ALLOCATABLE            :: IndexInactive(:)
    Vec                             :: TmpTplgy
    PetscBool                       :: Verbose
    CHARACTER(LEN=200)              :: ViewName,EventName

    CALL GetVerbose(Verbose,ierr)

    ! It obtains a temporal Fortran array where will be filled each 
    ! topology identifier.
    CALL DMCreateGlobalVector(DataMngr,TmpTplgy,ierr)
    CALL DMDAVecGetArrayF90(DataMngr,TmpTplgy,TmpTplgyArray,ierr)

    ! It quantifies Dirichlet, Neumman, and Cauchy boundary condition 
    ! on local and global processor.
    SizeTplgy(:)=0

    ! It fills the temporal Fortran array.

    ! It gets the global size from the geometry data manager.
    CALL DMDAGetInfo(DataMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
        & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,           &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,  &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,  &
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
                    SizeTplgy(2)=SizeTplgy(2)+1
                END IF
                TmpTplgyArray(i,j,k)=ValR
            END DO
        END DO
    END DO
    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,                  &
        & "[GetGeometry Event] WARNING: Topology File wasn't provi"//&
        & "ded. Default topology used is ...\n",ierr)

    ViewName="Tplgy"
    EventName="GetTopology"
    CALL ViewTopology(TmpTplgy,ViewName,EventName,ierr)

    ! It moves the temporal Fortran array to a petsc vector, Tplgy, 
    ! stored in Gmtry.
    CALL DMDAVecRestoreArrayF90(DataMngr,TmpTplgy,TmpTplgyArray,ierr)
    CALL DMCreateLocalVector(DataMngr,Tplgy,ierr)

    CALL DMGlobalToLocalBegin(DataMngr,TmpTplgy,INSERT_VALUES,Tplgy, &
        & ierr)
    CALL DMGlobalToLocalEnd(DataMngr,TmpTplgy,INSERT_VALUES,Tplgy,   &
        & ierr)

    CALL VecDestroy(TmpTplgy,ierr)

    CALL MPI_Bcast(SizeTplgy,4,MPI_INT,0,Comm,ierr)
    CALL MPI_Bcast(SizeInactive,1,MPI_INT,0,Comm,ierr)
    SizeTplgy(1)=widthG(1)*widthG(2)*widthG(3)

    ALLOCATE(IndexInactive(SizeInactive))
    CALL MPI_Bcast(IndexInactive,SizeInactive,MPI_INT,0,Comm,ierr)
    CALL ISCreateGeneral(Comm,SizeInactive,IndexInactive,            &
        & PETSC_COPY_VALUES,InactiveIS,ierr)
    DEALLOCATE(IndexInactive)

    IF (Verbose) CALL PetscSynchronizedPrintf(Comm,                  &
        & "[GetGeometry Event] Topology Identifiers were satisfact"//&
        & "orily created\n",ierr)

END SUBROUTINE GetTopology_3

SUBROUTINE UpdateTplgy(Gmtry,DirichIS,SourceIS,CauchyIS,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    TYPE(Geometry),INTENT(INOUT)        :: Gmtry
    IS,INTENT(IN)                       :: DirichIS,SourceIS,CauchyIS
    PetscErrorCode,INTENT(INOUT)        :: ierr

    PetscInt                            :: SizeDirich,SizeSource,    &
                                         & SizeCauchy
    PetscReal                           :: rTwo=2.d0,rThree=3.d0,    &
                                         & rFour=4.d0
    Vec                                 :: TmpTplgy,vTwo,vThree,vFour
    IS                                  :: NatOrderDirichIS,         &
                                         & NatOrderSourceIS,         &
                                         & NatOrderCauchyIS
    VecScatter                          :: DirichScatter,            &
                                         & SourceScatter,CauchyScatter


    ! Copying topology from the permanent Topology
    CALL VecCopy(Gmtry%pTplgy,Gmtry%Tplgy,ierr)

    CALL DMCreateGlobalVector(Gmtry%DataMngr,TmpTplgy,ierr)
    CALL DMLocalToGlobalBegin(Gmtry%DataMngr,Gmtry%Tplgy,            &
        & INSERT_VALUES,TmpTplgy,ierr)
    CALL DMLocalToGlobalEnd(Gmtry%DataMngr,Gmtry%Tplgy,INSERT_VALUES,&
        & TmpTplgy,ierr)

    CALL ISGetLocalSize(DirichIS,SizeDirich,ierr)
    CALL VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,SizeDirich,vTwo, &
        & ierr)
    CALL VecSet(vTwo,rTwo,ierr)
    CALL ISCreateStride(PETSC_COMM_WORLD,SizeDirich,0,1,             &
        & NatOrderDirichIS,ierr)
    CALL VecScatterCreate(vTwo,NatOrderDirichIS,TmpTplgy,DirichIS,   &
        & DirichScatter,ierr)
    CALL VecScatterBegin(DirichScatter,vTwo,TmpTplgy,INSERT_VALUES,  &
        & SCATTER_FORWARD,ierr)
    CALL VecScatterEnd(DirichScatter,vTwo,TmpTplgy,INSERT_VALUES,    &
        & SCATTER_FORWARD,ierr)
    CALL VecScatterDestroy(DirichScatter,ierr)
    CALL VecDestroy(vTwo,ierr)
    CALL ISDestroy(NatOrderDirichIS,ierr)

    CALL ISGetLocalSize(SourceIS,SizeSource,ierr)
    CALL VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,SizeSource,      & 
        & vThree,ierr)
    CALL VecSet(vThree,rThree,ierr)
    CALL ISCreateStride(PETSC_COMM_WORLD,SizeSource,0,1,             &
        & NatOrderSourceIS,ierr)
    CALL VecScatterCreate(vThree,NatOrderSourceIS,TmpTplgy,SourceIS, &
        & SourceScatter,ierr)
    CALL VecScatterBegin(SourceScatter,vThree,TmpTplgy,INSERT_VALUES,&
        & SCATTER_FORWARD,ierr)
    CALL VecScatterEnd(SourceScatter,vThree,TmpTplgy,INSERT_VALUES,  &
        & SCATTER_FORWARD,ierr)
    CALL VecScatterDestroy(SourceScatter,ierr)
    CALL VecDestroy(vThree,ierr)
    CALL ISDestroy(NatOrderSourceIS,ierr)

    CALL ISGetLocalSize(CauchyIS,SizeCauchy,ierr)
    CALL VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,SizeCauchy,vFour,&
        & ierr)
    CALL VecSet(vFour,rFour,ierr)
    CALL ISCreateStride(PETSC_COMM_WORLD,SizeCauchy,0,1,             &
        & NatOrderCauchyIS,ierr)
    CALL VecScatterCreate(vFour,NatOrderCauchyIS,TmpTplgy,CauchyIS,  &
        & CauchyScatter,ierr)
    CALL VecScatterBegin(CauchyScatter,vFour,TmpTplgy,INSERT_VALUES, &
        & SCATTER_FORWARD,ierr)
    CALL VecScatterEnd(CauchyScatter,vFour,TmpTplgy,INSERT_VALUES,   &
        & SCATTER_FORWARD,ierr)
    CALL VecScatterDestroy(CauchyScatter,ierr)
    CALL VecDestroy(vFour,ierr)
    CALL ISDestroy(NatOrderCauchyIS,ierr)

    CALL DMGlobalToLocalBegin(Gmtry%DataMngr,TmpTplgy,INSERT_VALUES, &
        & Gmtry%Tplgy,ierr)
    CALL DMGlobalToLocalEnd(Gmtry%DataMngr,TmpTplgy,INSERT_VALUES,   &
        & Gmtry%Tplgy,ierr)

END SUBROUTINE UpdateTplgy

SUBROUTINE GetLocalTopology(Gmtry,Ppt,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,Property

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

    PetscReal,POINTER                   :: TmpTplgyArray3D(:,:,:),   &
                                         & TmpTplgyArray2D(:,:),     &
                                         & xArray(:),yArray(:),      &
                                         & zArray(:)
    PetscInt                            :: i,j,k,widthG(3)
    PetscReal                           :: dx,dxB,dxF,dy,dyB,dyF,dz, &
                                         & dzB,dzF
    DMDAStencilType                     :: Stencil

    i=Ppt%Pstn%i
    j=Ppt%Pstn%j
    k=Ppt%Pstn%k

    ! It gets the global size from the geometry data manager.
    CALL DMDAGetInfo(Gmtry%DataMngr,PETSC_NULL_INTEGER,widthG(1),    & 
        & widthG(2),widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,  &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,  &
        & Stencil,ierr)

    dxB=0.D0;dx=0.D0;dxF=0.D0
    dyB=0.D0;dy=0.D0;dyF=0.D0
    dzB=0.D0;dz=0.D0;dzF=0.D0
    CALL VecGetArrayReadF90(Gmtry%x,xArray,ierr)
    IF (i.GT.0)             dxB=ABS(xArray(i+1)-xArray(i))
                            dx =ABS(xArray(i+2)-xArray(i+1))
    IF (i.LT.(widthG(1)-1)) dxF=ABS(xArray(i+3)-xArray(i+2))
    CALL VecRestoreArrayReadF90(Gmtry%x,xArray,ierr)

    CALL VecGetArrayReadF90(Gmtry%y,yArray,ierr)
    IF (j.GT.0)             dyB=ABS(yArray(j+1)-yArray(j))
                            dy =ABS(yArray(j+2)-yArray(j+1))
    IF (j.LT.(widthG(2)-1)) dyF=ABS(yArray(j+3)-yArray(j+2))
    CALL VecRestoreArrayReadF90(Gmtry%y,yArray,ierr)

    CALL VecGetArrayReadF90(Gmtry%z,zArray,ierr)
    IF (k.GT.0)             dzB=ABS(zArray(k+1)-zArray(k))
                            dz =ABS(zArray(k+2)-zArray(k+1))
    IF (k.LT.(widthG(3)-1)) dzF=ABS(zArray(k+3)-zArray(k+2))
    CALL VecRestoreArrayReadF90(Gmtry%z,zArray,ierr)

    ! Revisar los diferenciales de x,y,z que se usan como dividendo, 
    ! si los valores son cercanos a cero puede haber problemas!!!!

    IF (Stencil.EQ.DMDA_STENCIL_STAR) THEN
        Ppt%StnclType=1
    ELSEIF (Stencil.EQ.DMDA_STENCIL_BOX) THEN
        Ppt%StnclType=2
    END IF

    IF (widthG(3).NE.1) THEN
        CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,Gmtry%Tplgy,      &
            & TmpTplgyArray3D,ierr)
        IF (Ppt%StnclType.EQ.1) THEN
            ALLOCATE(Ppt%StnclTplgy(7),Ppt%dx(7),Ppt%dy(7),Ppt%dz(7))
            Ppt%StnclTplgy(:)=0
            Ppt%StnclTplgy(1)=INT(TmpTplgyArray3D(i  ,j  ,k-1))
            Ppt%StnclTplgy(2)=INT(TmpTplgyArray3D(i  ,j-1,k  ))
            Ppt%StnclTplgy(3)=INT(TmpTplgyArray3D(i-1,j  ,k  ))
            Ppt%StnclTplgy(4)=INT(TmpTplgyArray3D(i  ,j  ,k  ))
            Ppt%StnclTplgy(5)=INT(TmpTplgyArray3D(i+1,j  ,k  ))
            Ppt%StnclTplgy(6)=INT(TmpTplgyArray3D(i  ,j+1,k  ))
            Ppt%StnclTplgy(7)=INT(TmpTplgyArray3D(i  ,j  ,k+1))
            Ppt%dx(1)=dx   ; Ppt%dy(1)=dy   ; Ppt%dz(1)=dzB
            Ppt%dx(2)=dx   ; Ppt%dy(2)=dyB  ; Ppt%dz(2)=dz
            Ppt%dx(3)=dxB  ; Ppt%dy(3)=dy   ; Ppt%dz(3)=dz
            Ppt%dx(4)=dx   ; Ppt%dy(4)=dy   ; Ppt%dz(4)=dz
            Ppt%dx(5)=dxF  ; Ppt%dy(5)=dy   ; Ppt%dz(5)=dz
            Ppt%dx(6)=dx   ; Ppt%dy(6)=dyF  ; Ppt%dz(6)=dz
            Ppt%dx(7)=dx   ; Ppt%dy(7)=dy   ; Ppt%dz(7)=dzF
        ELSEIF (Ppt%StnclType.EQ.2) THEN
            ALLOCATE(Ppt%StnclTplgy(19),Ppt%dx(19),Ppt%dy(19),       &
                & Ppt%dz(19))
            Ppt%StnclTplgy(:)=0
            Ppt%StnclTplgy(1)= INT(TmpTplgyArray3D(i  ,j-1,k-1))
            Ppt%StnclTplgy(2)= INT(TmpTplgyArray3D(i-1,j  ,k-1))
            Ppt%StnclTplgy(3)= INT(TmpTplgyArray3D(i  ,j  ,k-1))
            Ppt%StnclTplgy(4)= INT(TmpTplgyArray3D(i+1,j  ,k-1))
            Ppt%StnclTplgy(5)= INT(TmpTplgyArray3D(i  ,j+1,k-1))
            Ppt%StnclTplgy(6)= INT(TmpTplgyArray3D(i-1,j-1,k  ))
            Ppt%StnclTplgy(7)= INT(TmpTplgyArray3D(i  ,j-1,k  ))
            Ppt%StnclTplgy(8)= INT(TmpTplgyArray3D(i+1,j-1,k  ))
            Ppt%StnclTplgy(9)= INT(TmpTplgyArray3D(i-1,j  ,k  ))
            Ppt%StnclTplgy(10)=INT(TmpTplgyArray3D(i  ,j  ,k  ))
            Ppt%StnclTplgy(11)=INT(TmpTplgyArray3D(i+1,j  ,k  ))
            Ppt%StnclTplgy(12)=INT(TmpTplgyArray3D(i-1,j+1,k  ))
            Ppt%StnclTplgy(13)=INT(TmpTplgyArray3D(i  ,j+1,k  ))
            Ppt%StnclTplgy(14)=INT(TmpTplgyArray3D(i+1,j+1,k  ))
            Ppt%StnclTplgy(15)=INT(TmpTplgyArray3D(i  ,j-1,k+1))
            Ppt%StnclTplgy(16)=INT(TmpTplgyArray3D(i-1,j  ,k+1))
            Ppt%StnclTplgy(17)=INT(TmpTplgyArray3D(i  ,j  ,k+1))
            Ppt%StnclTplgy(18)=INT(TmpTplgyArray3D(i+1,j  ,k+1))
            Ppt%StnclTplgy(19)=INT(TmpTplgyArray3D(i  ,j+1,k+1))
            Ppt%dx(1)= dx  ; Ppt%dy(1)= dyB ; Ppt%dz(1)= dzB
            Ppt%dx(2)= dxB ; Ppt%dy(2)= dy  ; Ppt%dz(2)= dzB
            Ppt%dx(3)= dx  ; Ppt%dy(3)= dy  ; Ppt%dz(3)= dzB
            Ppt%dx(4)= dxF ; Ppt%dy(4)= dy  ; Ppt%dz(4)= dzB
            Ppt%dx(5)= dx  ; Ppt%dy(5)= dyF ; Ppt%dz(5)= dzB
            Ppt%dx(6)= dxB ; Ppt%dy(6)= dyB ; Ppt%dz(6)= dz
            Ppt%dx(7)= dx  ; Ppt%dy(7)= dyB ; Ppt%dz(7)= dz
            Ppt%dx(8)= dxF ; Ppt%dy(8)= dyB ; Ppt%dz(8)= dz
            Ppt%dx(9)= dxB ; Ppt%dy(9)= dy  ; Ppt%dz(9)= dz
            Ppt%dx(10)=dx  ; Ppt%dy(10)=dy  ; Ppt%dz(10)=dz
            Ppt%dx(11)=dxF ; Ppt%dy(11)=dy  ; Ppt%dz(11)=dz
            Ppt%dx(12)=dxB ; Ppt%dy(12)=dyF ; Ppt%dz(12)=dz
            Ppt%dx(13)=dx  ; Ppt%dy(13)=dyF ; Ppt%dz(13)=dz
            Ppt%dx(14)=dxF ; Ppt%dy(14)=dyF ; Ppt%dz(14)=dz
            Ppt%dx(15)=dx  ; Ppt%dy(15)=dyB ; Ppt%dz(15)=dzF
            Ppt%dx(16)=dxB ; Ppt%dy(16)=dy  ; Ppt%dz(16)=dzF
            Ppt%dx(17)=dx  ; Ppt%dy(17)=dy  ; Ppt%dz(17)=dzF
            Ppt%dx(18)=dxF ; Ppt%dy(18)=dy  ; Ppt%dz(18)=dzF
            Ppt%dx(19)=dx  ; Ppt%dy(19)=dyF ; Ppt%dz(19)=dzF
        ELSE
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,           &
                & "[ERROR] Ppt%StnclType must be 1 or 2.\n",ierr)
            STOP
        END IF
    ELSE
        CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,Gmtry%Tplgy,      &
            & TmpTplgyArray2D,ierr)
        IF (Ppt%StnclType.EQ.1) THEN
            ALLOCATE(Ppt%StnclTplgy(7),Ppt%dx(7),Ppt%dy(7),Ppt%dz(7))
            Ppt%StnclTplgy(:)=0
            Ppt%StnclTplgy(1)=INT(TmpTplgyArray2D(i  ,j  ))
            Ppt%StnclTplgy(2)=INT(TmpTplgyArray2D(i  ,j-1))
            Ppt%StnclTplgy(3)=INT(TmpTplgyArray2D(i-1,j  ))
            Ppt%StnclTplgy(4)=INT(TmpTplgyArray2D(i  ,j  ))
            Ppt%StnclTplgy(5)=INT(TmpTplgyArray2D(i+1,j  ))
            Ppt%StnclTplgy(6)=INT(TmpTplgyArray2D(i  ,j+1))
            Ppt%StnclTplgy(7)=INT(TmpTplgyArray2D(i  ,j  ))
            Ppt%dx(1)=dx   ; Ppt%dy(1)=dy   ; Ppt%dz(1)=dzB
            Ppt%dx(2)=dx   ; Ppt%dy(2)=dyB  ; Ppt%dz(2)=dz
            Ppt%dx(3)=dxB  ; Ppt%dy(3)=dy   ; Ppt%dz(3)=dz
            Ppt%dx(4)=dx   ; Ppt%dy(4)=dy   ; Ppt%dz(4)=dz
            Ppt%dx(5)=dxF  ; Ppt%dy(5)=dy   ; Ppt%dz(5)=dz
            Ppt%dx(6)=dx   ; Ppt%dy(6)=dyF  ; Ppt%dz(6)=dz
            Ppt%dx(7)=dx   ; Ppt%dy(7)=dy   ; Ppt%dz(7)=dzF
        ELSEIF (Ppt%StnclType.EQ.2) THEN
            ALLOCATE(Ppt%StnclTplgy(19),Ppt%dx(19),Ppt%dy(19),       &
                & Ppt%dz(19))
            Ppt%StnclTplgy(:)=0
            Ppt%StnclTplgy(1)= INT(TmpTplgyArray2D(i  ,j-1))
            Ppt%StnclTplgy(2)= INT(TmpTplgyArray2D(i-1,j  ))
            Ppt%StnclTplgy(3)= INT(TmpTplgyArray2D(i  ,j  ))
            Ppt%StnclTplgy(4)= INT(TmpTplgyArray2D(i+1,j  ))
            Ppt%StnclTplgy(5)= INT(TmpTplgyArray2D(i  ,j+1))
            Ppt%StnclTplgy(6)= INT(TmpTplgyArray2D(i-1,j-1))
            Ppt%StnclTplgy(7)= INT(TmpTplgyArray2D(i  ,j-1))
            Ppt%StnclTplgy(8)= INT(TmpTplgyArray2D(i+1,j-1))
            Ppt%StnclTplgy(9)= INT(TmpTplgyArray2D(i-1,j  ))
            Ppt%StnclTplgy(10)=INT(TmpTplgyArray2D(i  ,j  ))
            Ppt%StnclTplgy(11)=INT(TmpTplgyArray2D(i+1,j  ))
            Ppt%StnclTplgy(12)=INT(TmpTplgyArray2D(i-1,j+1))
            Ppt%StnclTplgy(13)=INT(TmpTplgyArray2D(i  ,j+1))
            Ppt%StnclTplgy(14)=INT(TmpTplgyArray2D(i+1,j+1))
            Ppt%StnclTplgy(15)=INT(TmpTplgyArray2D(i  ,j-1))
            Ppt%StnclTplgy(16)=INT(TmpTplgyArray2D(i-1,j  ))
            Ppt%StnclTplgy(17)=INT(TmpTplgyArray2D(i  ,j  ))
            Ppt%StnclTplgy(18)=INT(TmpTplgyArray2D(i+1,j  ))
            Ppt%StnclTplgy(19)=INT(TmpTplgyArray2D(i  ,j+1))
            Ppt%dx(1)= dx  ; Ppt%dy(1)= dyB ; Ppt%dz(1)= dzB
            Ppt%dx(2)= dxB ; Ppt%dy(2)= dy  ; Ppt%dz(2)= dzB
            Ppt%dx(3)= dx  ; Ppt%dy(3)= dy  ; Ppt%dz(3)= dzB
            Ppt%dx(4)= dxF ; Ppt%dy(4)= dy  ; Ppt%dz(4)= dzB
            Ppt%dx(5)= dx  ; Ppt%dy(5)= dyF ; Ppt%dz(5)= dzB
            Ppt%dx(6)= dxB ; Ppt%dy(6)= dyB ; Ppt%dz(6)= dz
            Ppt%dx(7)= dx  ; Ppt%dy(7)= dyB ; Ppt%dz(7)= dz
            Ppt%dx(8)= dxF ; Ppt%dy(8)= dyB ; Ppt%dz(8)= dz
            Ppt%dx(9)= dxB ; Ppt%dy(9)= dy  ; Ppt%dz(9)= dz
            Ppt%dx(10)=dx  ; Ppt%dy(10)=dy  ; Ppt%dz(10)=dz
            Ppt%dx(11)=dxF ; Ppt%dy(11)=dy  ; Ppt%dz(11)=dz
            Ppt%dx(12)=dxB ; Ppt%dy(12)=dyF ; Ppt%dz(12)=dz
            Ppt%dx(13)=dx  ; Ppt%dy(13)=dyF ; Ppt%dz(13)=dz
            Ppt%dx(14)=dxF ; Ppt%dy(14)=dyF ; Ppt%dz(14)=dz
            Ppt%dx(15)=dx  ; Ppt%dy(15)=dyB ; Ppt%dz(15)=dzF
            Ppt%dx(16)=dxB ; Ppt%dy(16)=dy  ; Ppt%dz(16)=dzF
            Ppt%dx(17)=dx  ; Ppt%dy(17)=dy  ; Ppt%dz(17)=dzF
            Ppt%dx(18)=dxF ; Ppt%dy(18)=dy  ; Ppt%dz(18)=dzF
            Ppt%dx(19)=dx  ; Ppt%dy(19)=dyF ; Ppt%dz(19)=dzF
        ELSE
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,           &
                & "[ERROR] Ppt%StnclType must be 1 or 2.\n",ierr)
            STOP
        END IF
        CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,Gmtry%Tplgy,  &
            & TmpTplgyArray2D,ierr)
    END IF

END SUBROUTINE GetLocalTopology

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
    ! Changing temporal vector from application ordering to PETSc 
    ! ordering
    CALL DMDAGetAO(DataMngr,AppOrd,ierr)
    CALL VecGetOwnershipRange(PetscVec,RangeLow,RangeHigh,ierr)
    CALL ISCreateGeneral(PETSC_COMM_WORLD,RangeHigh-RangeLow,        &
        & (/(i,i=RangeLow,RangeHigh)/),PETSC_COPY_VALUES,PetscIS,ierr)
    CALL ISDuplicate(PetscIS,AppIS,ierr)
    CALL ISCopy(PetscIS,AppIS,ierr)
    CALL AOPetscToApplicationIS(AppOrd,AppIS,ierr)
    CALL VecScatterCreate(AppVec,AppIS,PetscVec,PetscIS,Scatter,ierr)
    CALL VecScatterBegin(Scatter,AppVec,PetscVec,INSERT_VALUES,      & 
        & SCATTER_FORWARD,ierr)
    CALL VecScatterEnd(Scatter,AppVec,PetscVec,INSERT_VALUES,        &
        & SCATTER_FORWARD,ierr)

    CALL VecScatterDestroy(Scatter,ierr)
    CALL ISDestroy(AppIS,ierr)
    CALL ISDestroy(PetscIS,ierr)
    CALL VecCopy(PetscVec,AppVec,ierr)
    CALL VecDestroy(PetscVec,ierr)

END SUBROUTINE VecApplicationToPetsc

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
    PetscLogDouble                  :: EventFlops=0.d0

    ClassName="Geometry"
    CALL PetscClassIdRegister(ClassName,ClassID,ierr)
    EventName="DestroyGeometry"
    CALL PetscLogEventRegister(EventName,ClassID,Event,ierr)
    CALL PetscLogEventBegin(Event,PETSC_NULL_OBJECT,                 &
        & PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,      &
        & "["//ADJUSTL(TRIM(EventName))//" Event] Inizialited\n",ierr)

    CALL DMDestroy(Gmtry%DataMngr,ierr)
    CALL VecDestroy(Gmtry%Tplgy,ierr)
    CALL VecDestroy(Gmtry%x,ierr)
    CALL VecDestroy(Gmtry%y,ierr)
    CALL VecDestroy(Gmtry%z,ierr)

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,      &
        & "["//ADJUSTL(TRIM(EventName))//" Event] Finalized\n",ierr)
    
    CALL PetscLogFlops(EventFlops,ierr)
    CALL PetscLogEventEnd(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
        & PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

END SUBROUTINE DestroyGeometry

END MODULE ANISOFLOW_Geometry