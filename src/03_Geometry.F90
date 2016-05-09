MODULE ANISOFLOW_Geometry

! ANISOFLOW_Geometry it's a module that contains routines to manage geometry variables.

    USE ANISOFLOW_Types, ONLY : Geometry

    IMPLICIT NONE

CONTAINS

 !  - GetGeometry: It's a routine that fills a Geometry data structure with input files provided by
 !                 the user.
 !    > OUT: Gmtry, ierr.
 !      + Gmtry: It's a Geometry data structure filled with input files provided by the user.
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.
 !    > NOTES: The Geometry is filled in two stages: to create a Data Manager (DataMngr) to control
 !             information related to a regular rectangular grid, and Topology that describes the
 !             geometry with identifiers.

SUBROUTINE GetGeometry(Gmtry,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    TYPE(Geometry),INTENT(OUT)      :: Gmtry

    CALL GetDataMngr(Gmtry%DataMngr,ierr)
    CALL GetGridCoord(Gmtry%DataMngr,Gmtry%x,Gmtry%y,Gmtry%z,ierr)
    CALL GetTopology(Gmtry%DataMngr,Gmtry%Tplgy,Gmtry%DirichIS,Gmtry%NeummanIS,Gmtry%CauchyIS,ierr)

END SUBROUTINE GetGeometry

 !  - GetDataMngr: It's a routine that creates and fills the information related with a regular 
 !                 rectangular grid.
 !    > OUT: DataMngr, ierr.
 !      + DataMngr: It's a DMDA PETSc structure that stores the information related with a regular 
 !                 rectangular grid.
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.
 !    > NOTES: The DataMngr takes the size of the domain and assigning a subdomain to each processor.

SUBROUTINE GetDataMngr(DataMngr,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    DM,INTENT(INOUT)                :: DataMngr

    PetscMPIInt                     :: process
    PetscInt                        :: u,widthG(3)
    CHARACTER(LEN=200)              :: InputDir,InputFileGmtry,Route
    TYPE(InputTypeVar)              :: InputType
    TYPE(RunOptionsVar)             :: RunOptions
    DMDAStencilType                 :: Stencil

    PARAMETER(u=01)

    ! It obtains the route to open a geometry file.
    CALL GetInputDir(InputDir,ierr)
    CALL GetInputType(InputType,ierr)
    CALL GetInputFileGmtry(InputFileGmtry,ierr)

    ! It obtains run options.
    CALL GetRunOptions(RunOptions,ierr)

    CALL MPI_Comm_rank(MPI_COMM_WORLD,process,ierr)

    ! It obtains the global size of the domain on the first processor.
    IF (process.EQ.0) THEN
        Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFileGmtry))
        OPEN(u,FILE=TRIM(Route),STATUS='OLD',ACTION='READ')

        ! It gets the size of the domain depending on the geometry file input.
        !   1: Defined by Blessent. An example is provided in "../ex/Blessent/in/tsim_USMH.asc"
        !   2: Defined by Perez. An example is provided in "../ex/Perez/in/sanpck.domnRST"
        IF (InputType%Gmtry.EQ.1) THEN
            READ(u, '((I10),(I10),(I10))')widthG(1),widthG(2),widthG(3)
        END IF
        CLOSE(u)
    END IF

    ! It broadcasts the global size to other processors.
    CALL MPI_Bcast(widthG,3,MPI_INT, 0, PETSC_COMM_WORLD,ierr)

    ! It decides the stencil shape depending on the scheme used.
    IF (RunOptions%Scheme.EQ.1) THEN
        Stencil=DMDA_STENCIL_STAR
    ELSEIF (RunOptions%Scheme.EQ.2) THEN
        Stencil=DMDA_STENCIL_BOX
    ELSE
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "ERROR: Run_options_scheme command must be an integer between 1 and 2\n",ierr)
        STOP
    END IF

    ! It creates the Data Manager to Distributed Arrays by the information provided.
    CALL DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,      &
        & DM_BOUNDARY_GHOSTED,Stencil,widthG(1),widthG(2),widthG(3),     &
        & PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,PETSC_NULL_INTEGER,       &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,DataMngr,ierr)

END SUBROUTINE GetDataMngr

SUBROUTINE GetGridCoord(DataMngr,x,y,z,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    DM,INTENT(IN)                   :: DataMngr
    Vec,INTENT(OUT)                 :: x,y,z

    CHARACTER(LEN=200)              :: InputDir,InputFileGmtry,Route
    TYPE(InputTypeVar)              :: InputType
    PetscInt                        :: widthG(3),size,i
    PetscReal                       :: Value


    ! It obtains the route to open a geometry file.
    CALL GetInputDir(InputDir,ierr)
    CALL GetInputType(InputType,ierr)
    CALL GetInputFileGmtry(InputFileGmtry,ierr)

    ! It gets the global size from the geometry data manager.
    CALL DMDAGetInfo(DataMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
        & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,ierr)

    size=widthG(1)+1
    CALL VecCreateSeq(PETSC_COMM_SELF,size,x,ierr)

    size=widthG(2)+1
    CALL VecCreateSeq(PETSC_COMM_SELF,size,y,ierr)

    size=widthG(3)+1
    CALL VecCreateSeq(PETSC_COMM_SELF,size,z,ierr)

    IF (InputType%Gmtry.EQ.1) THEN
        ! Default DX=DY=DZ=1.0
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
    END IF

    CALL VecAssemblyBegin(x,ierr)
    CALL VecAssemblyEnd(x,ierr)

    CALL VecAssemblyBegin(y,ierr)
    CALL VecAssemblyEnd(y,ierr)
    
    CALL VecAssemblyBegin(z,ierr)
    CALL VecAssemblyEnd(z,ierr)

END SUBROUTINE GetGridCoord

 !  - GetTopology: It's a routine that creates and fills the information related to topology.
 !                 It creates a vector and index sets to describe the geometry. 
 !    > OUT: Gmtry, ierr.
 !      + Gmtry: It's a Geometry data structure filled with input files provided by the user.
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.


SUBROUTINE GetTopology(DataMngr,Tplgy,DirichIS,NeummanIS,CauchyIS,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscdmda.h90>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    DM,INTENT(IN)                   :: DataMngr
    Vec,INTENT(OUT)                 :: Tplgy
    IS,INTENT(OUT)                  :: DirichIS,NeummanIS,CauchyIS

    PetscReal,POINTER               :: TmpTplgyArray(:,:,:)
    PetscReal                       :: ValR
    PetscInt                        :: i,j,k,widthL(3),widthG(3),corn(3),BCLenL(3),BCLenG(3)
    CHARACTER(LEN=200)              :: InputDir,InputFileGmtry,Route
    TYPE(InputTypeVar)              :: InputType
    Vec                             :: TmpTplgy

    ! It obtains the route to open a geometry file.
    CALL GetInputDir(InputDir,ierr)
    CALL GetInputType(InputType,ierr)
    CALL GetInputFileGmtry(InputFileGmtry,ierr)

    ! It obtains a temporal Fortran array where will be filled each topology identifier.
    CALL DMCreateGlobalVector(DataMngr,TmpTplgy,ierr)
    CALL DMDAVecGetArrayF90(DataMngr,TmpTplgy,TmpTplgyArray,ierr)

    ! It quantifies Dirichlet, Neumman, and Cauchy boundary condition on local and global processor.
    BCLenL(:)=0
    BCLenG(:)=0

    ! It fills the temporal Fortran array depending on the topology file input.
    !   1: Default topology, it doesn't need a file. The first border layer is Dirichlet, active in other case
    IF (InputType%Tplgy.EQ.1) THEN

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
    END IF

    ! It moves the temporal Fortran array to a petsc vector, Tplgy, stored in Gmtry.
    CALL DMDAVecRestoreArrayF90(DataMngr,TmpTplgy,TmpTplgyArray,ierr)
    CALL DMCreateLocalVector(DataMngr,Tplgy,ierr)

    CALL DMGlobalToLocalBegin(DataMngr,TmpTplgy,INSERT_VALUES,  &
        & Tplgy,ierr)
    CALL DMGlobalToLocalEnd(DataMngr,TmpTplgy,INSERT_VALUES,    &
        & Tplgy,ierr)

    CALL VecDestroy(TmpTplgy,ierr)

    CALL MPI_ALLREDUCE(BCLenL,BCLenG,3,MPI_INT,MPI_SUM,PETSC_COMM_WORLD,ierr)

    CALL GetTopologyBC(DataMngr,Tplgy,BCLenG,DirichIS,NeummanIS,CauchyIS,ierr)

END SUBROUTINE GetTopology

SUBROUTINE GetTopologyBC(DataMngr,Tplgy,BCLenG,DirichIS,NeummanIS,CauchyIS,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscis.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    DM,INTENT(IN)                   :: DataMngr
    Vec,INTENT(IN)                  :: Tplgy
    PetscInt,INTENT(IN)             :: BCLenG(3)
    IS,INTENT(OUT)                  :: DirichIS,NeummanIS,CauchyIS

    PetscInt,ALLOCATABLE            :: IndexDirich(:),IndexNeumman(:)
    PetscInt,ALLOCATABLE            :: IndexCauchy(:)
    PetscInt                        :: i,Low,High,Size,ValI,Count(3),SumaL(3),SumaG(3)
    PetscReal,POINTER               :: TmpTplgyArray(:)
    Vec                             :: TmpTplgy

    ALLOCATE(IndexDirich(BCLenG(1)))
    ALLOCATE(IndexNeumman(BCLenG(2)))
    ALLOCATE(IndexCauchy(BCLenG(3)))
    IndexDirich(:)=0
    IndexNeumman(:)=0
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
            ELSEIF ((ValI.EQ.3).OR.(ValI.EQ.4).OR.(ValI.EQ.5)) THEN
                IndexNeumman(Count(2))=i
                SumaL(2)=1
            ELSEIF (ValI.EQ.6) THEN
                IndexCauchy(Count(3))=i
                SumaL(3)=1
            END IF
        END IF
        CALL MPI_ALLREDUCE(SumaL,SumaG,3,MPI_INT,MPI_SUM,PETSC_COMM_WORLD,ierr)
        Count(:)=Count(:)+SumaG(:)
    END DO

    CALL VecRestoreArrayReadF90(TmpTplgy,TmpTplgyArray,ierr)

    CALL DMGlobalToLocalBegin(DataMngr,TmpTplgy,INSERT_VALUES,  &
        & Tplgy,ierr)
    CALL DMGlobalToLocalEnd(DataMngr,TmpTplgy,INSERT_VALUES,    &
        & Tplgy,ierr)

    CALL VecDestroy(TmpTplgy,ierr)

    CALL ISCreateGeneral(MPI_COMM_WORLD,BCLenG(1),IndexDirich,                 &
        & PETSC_COPY_VALUES,DirichIS,ierr)
    CALL ISCreateGeneral(PETSC_COMM_WORLD,BCLenG(2),IndexNeumman,              &
        & PETSC_COPY_VALUES,NeummanIS,ierr)
    CALL ISCreateGeneral(PETSC_COMM_WORLD,BCLenG(3),IndexCauchy,               &
        & PETSC_COPY_VALUES,CauchyIS,ierr)

    DEALLOCATE(IndexDirich)
    DEALLOCATE(IndexNeumman)
    DEALLOCATE(IndexCauchy)

END SUBROUTINE GetTopologyBC

SUBROUTINE GetLocalTopology(Gmtry,Ppt,ierr)

    USE ANISOFLOW_Types, ONLY : Property
    USE ANISOFLOW_Interface

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
    Ppt%dzB=zArray(k+1)-zArray(k)
    Ppt%dz =zArray(k+2)-zArray(k+1)
    Ppt%dzF=zArray(k+3)-zArray(k+2)
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
            & "ERROR: Run_options_scheme command must be an integer between 1 and 2\n",ierr)
        STOP
    END IF
    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,Gmtry%Tplgy,TmpTplgyArray,ierr)

END SUBROUTINE GetLocalTopology

SUBROUTINE GeometryDestroy(Gmtry,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    TYPE(Geometry),INTENT(IN)       :: Gmtry ! INOUT or IN?

    CALL DMDestroy(Gmtry%DataMngr,ierr)
    CALL VecDestroy(Gmtry%Tplgy,ierr)
    CALL VecDestroy(Gmtry%x,ierr)
    CALL VecDestroy(Gmtry%y,ierr)
    CALL VecDestroy(Gmtry%z,ierr)

END SUBROUTINE GeometryDestroy

END MODULE ANISOFLOW_Geometry