MODULE ANISOFLOW_Geometry

    USE ANISOFLOW_Types, ONLY : Geometry

    IMPLICIT NONE

CONTAINS

SUBROUTINE GetGeometry(Gmtry,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    TYPE(Geometry),INTENT(OUT)      :: Gmtry

    CALL GetDstrMngr(Gmtry,ierr)
    CALL GetTopology(Gmtry,ierr)

END SUBROUTINE GetGeometry


SUBROUTINE GetDstrMngr(Gmtry,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    TYPE(Geometry),INTENT(INOUT)    :: Gmtry

    PetscMPIInt                     :: process
    PetscInt                        :: u,widthG(3)
    CHARACTER(LEN=200)              :: InputDir,InputFileGmtry,Route
    TYPE(InputTypeVar)              :: InputType
    TYPE(RunOptionsVar)             :: RunOptions
    DMDAStencilType                 :: Stencil

    PARAMETER(u=01)

    CALL GetInputDir(InputDir,ierr)
    CALL GetInputType(InputType,ierr)
    CALL GetInputFileGmtry(InputFileGmtry,ierr)

    CALL GetRunOptions(RunOptions,ierr)

    ! Information about geometry is read from master processor and sent to slaves

    CALL MPI_Comm_rank(MPI_COMM_WORLD,process,ierr)

    IF (process.EQ.0) THEN
        Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFileGmtry))
        OPEN(u,FILE=TRIM(Route),STATUS='OLD',ACTION='READ')
        IF (InputType%Gmtry.EQ.1) THEN
            READ(u, '((I10),(I10),(I10))')widthG(1),widthG(2),widthG(3)
        END IF
        CLOSE(u)
    END IF
    CALL MPI_Bcast(widthG,3,MPI_INT, 0, PETSC_COMM_WORLD,ierr)

    IF (RunOptions%Scheme.EQ.1) THEN
        Stencil=DMDA_STENCIL_STAR
    ELSEIF (RunOptions%Scheme.EQ.2) THEN
        Stencil=DMDA_STENCIL_BOX
    ELSE
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "ERROR: Run_options_scheme command must be an integer between 1 and 2\n",ierr)
        STOP
    END IF

    CALL DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,      &
        & DM_BOUNDARY_GHOSTED,Stencil,widthG(1),widthG(2),widthG(3),     &
        & PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,PETSC_NULL_INTEGER,       &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,Gmtry%DstrMngr,ierr)

END SUBROUTINE GetDstrMngr

SUBROUTINE GetTopology(Gmtry,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscdmda.h90>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    TYPE(Geometry),INTENT(INOUT)    :: Gmtry

    PetscReal,POINTER               :: TmpTplgyArray(:,:,:)
    PetscReal                       :: ValR
    PetscInt                        :: i,j,k,widthL(3),widthG(3),corn(3),BCLenL(3),BCLenG(3)
    CHARACTER(LEN=200)              :: InputDir,InputFileGmtry,Route
    TYPE(InputTypeVar)              :: InputType
    Vec                             :: TmpTplgy

    CALL GetInputDir(InputDir,ierr)
    CALL GetInputType(InputType,ierr)
    CALL GetInputFileGmtry(InputFileGmtry,ierr)

    CALL DMCreateGlobalVector(Gmtry%DstrMngr,TmpTplgy,ierr)

    CALL DMDAVecGetArrayF90(Gmtry%DstrMngr,TmpTplgy,TmpTplgyArray,ierr)

    CALL DMDAGetCorners(Gmtry%DstrMngr,corn(1),corn(2),corn(3),widthL(1),      &
        & widthL(2),widthL(3),ierr)
    BCLenL(:)=0
    BCLenG(:)=0

    IF (InputType%Tplgy.EQ.1) THEN
        ! Default topology

        CALL DMDAGetInfo(Gmtry%DstrMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
            & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
            & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
            & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
            & PETSC_NULL_INTEGER,ierr)

        ValR=1.0 ! Active
        TmpTplgyArray(:,:,:)=ValR
        DO k=corn(3),corn(3)+widthL(3)-1
            DO j=corn(2),corn(2)+widthL(2)-1
                DO i=corn(1),corn(1)+widthL(1)-1
                    ValR=1.0 ! Active
                    IF ((k.EQ.0).AND.((i.EQ.0).OR.(i.EQ.(widthG(1)-1)).OR.     &
                    & (j.EQ.0).OR.(j.EQ.(widthG(2)-1)))) THEN
                        ValR=2.0    ! Dirichlet
                        BCLenL(1)=BCLenL(1)+1
                    ELSEIF ((i.EQ.0).OR.(i.EQ.(widthG(1)-1))) THEN
                        ValR=3.0    ! Neumann in x
                        BCLenL(2)=BCLenL(2)+1
                    ELSEIF ((j.EQ.0).OR.(j.EQ.(widthG(2)-1))) THEN
                        ValR=4.0    ! Neumann in y
                        BCLenL(2)=BCLenL(2)+1
                    ELSEIF ((k.EQ.0).OR.(k.EQ.(widthG(3)-1))) THEN
                        ValR=5.0    ! Neumann in z
                        BCLenL(2)=BCLenL(2)+1
                    END IF
                    TmpTplgyArray(i,j,k)=ValR
                END DO
            END DO
        END DO
    END IF

    CALL DMDAVecRestoreArrayF90(Gmtry%DstrMngr,TmpTplgy,TmpTplgyArray,ierr)

    CALL DMCreateLocalVector(Gmtry%DstrMngr,Gmtry%Tplgy,ierr)

    CALL DMGlobalToLocalBegin(Gmtry%DstrMngr,TmpTplgy,INSERT_VALUES,  &
        & Gmtry%Tplgy,ierr)
    CALL DMGlobalToLocalEnd(Gmtry%DstrMngr,TmpTplgy,INSERT_VALUES,    &
        & Gmtry%Tplgy,ierr)

    CALL VecDestroy(TmpTplgy,ierr)

    CALL MPI_ALLREDUCE(BCLenL,BCLenG,3,MPI_INT,MPI_SUM,PETSC_COMM_WORLD,ierr)
    CALL GetTopologyBC(Gmtry,BCLenG,ierr)

END SUBROUTINE GetTopology

SUBROUTINE GetTopologyBC(Gmtry,BCLenG,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscis.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    TYPE(Geometry),INTENT(INOUT)    :: Gmtry
    PetscInt,INTENT(IN)             :: BCLenG(3)

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
    CALL DMCreateGlobalVector(Gmtry%DstrMngr,TmpTplgy,ierr)

    CALL DMLocalToGlobalBegin(Gmtry%DstrMngr,Gmtry%Tplgy,INSERT_VALUES,  &
        & TmpTplgy,ierr)
    CALL DMLocalToGlobalEnd(Gmtry%DstrMngr,Gmtry%Tplgy,INSERT_VALUES,    &
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

    CALL DMGlobalToLocalBegin(Gmtry%DstrMngr,TmpTplgy,INSERT_VALUES,  &
        & Gmtry%Tplgy,ierr)
    CALL DMGlobalToLocalEnd(Gmtry%DstrMngr,TmpTplgy,INSERT_VALUES,    &
        & Gmtry%Tplgy,ierr)

    CALL VecDestroy(TmpTplgy,ierr)

    CALL ISCreateGeneral(MPI_COMM_WORLD,BCLenG(1),IndexDirich,                 &
        & PETSC_COPY_VALUES,Gmtry%DirichIS,ierr)
    CALL ISCreateGeneral(PETSC_COMM_WORLD,BCLenG(2),IndexNeumman,              &
        & PETSC_COPY_VALUES,Gmtry%NeummanIS,ierr)
    CALL ISCreateGeneral(PETSC_COMM_WORLD,BCLenG(3),IndexCauchy,               &
        & PETSC_COPY_VALUES,Gmtry%CauchyIS,ierr)


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
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscdmda.h90>

    TYPE(Geometry),INTENT(IN)           :: Gmtry
    TYPE(Property),INTENT(INOUT)        :: Ppt
    PetscErrorCode,INTENT(INOUT)        :: ierr

    PetscReal,POINTER                   :: TmpTplgyArray(:,:,:)
    PetscInt                            :: i,j,k
    TYPE(RunOptionsVar)                 :: RunOptions

    i=Ppt%Pstn%i
    j=Ppt%Pstn%j
    k=Ppt%Pstn%k

    CALL GetRunOptions(RunOptions,ierr)

    CALL DMDAVecGetArrayReadF90(Gmtry%DstrMngr,Gmtry%Tplgy,TmpTplgyArray,ierr)

    IF (RunOptions%Scheme.EQ.1) THEN
        ALLOCATE(Ppt%Tplgy(7))
        Ppt%Tplgy(:)=0
        Ppt%Tplgy(1)=INT(TmpTplgyArray(i  ,j  ,k-1))
        Ppt%Tplgy(2)=INT(TmpTplgyArray(i  ,j-1,k  ))
        Ppt%Tplgy(3)=INT(TmpTplgyArray(i-1,j  ,k  ))
        Ppt%Tplgy(4)=INT(TmpTplgyArray(i  ,j  ,k  ))
        Ppt%Tplgy(5)=INT(TmpTplgyArray(i+1,j  ,k  ))
        Ppt%Tplgy(6)=INT(TmpTplgyArray(i  ,j+1,k  ))
        Ppt%Tplgy(7)=INT(TmpTplgyArray(i  ,j  ,k+1))
    ELSEIF (RunOptions%Scheme.EQ.2) THEN
        ALLOCATE(Ppt%Tplgy(19))
        Ppt%Tplgy(:)=0
        Ppt%Tplgy(1)= INT(TmpTplgyArray(i  ,j-1,k-1))
        Ppt%Tplgy(2)= INT(TmpTplgyArray(i-1,j  ,k-1))
        Ppt%Tplgy(3)= INT(TmpTplgyArray(i  ,j  ,k-1))
        Ppt%Tplgy(4)= INT(TmpTplgyArray(i+1,j  ,k-1))
        Ppt%Tplgy(5)= INT(TmpTplgyArray(i  ,j+1,k-1))
        Ppt%Tplgy(6)= INT(TmpTplgyArray(i-1,j-1,k  ))
        Ppt%Tplgy(7)= INT(TmpTplgyArray(i  ,j-1,k  ))
        Ppt%Tplgy(8)= INT(TmpTplgyArray(i+1,j-1,k  ))
        Ppt%Tplgy(9)= INT(TmpTplgyArray(i-1,j  ,k  ))
        Ppt%Tplgy(10)=INT(TmpTplgyArray(i  ,j  ,k  ))
        Ppt%Tplgy(11)=INT(TmpTplgyArray(i+1,j  ,k  ))
        Ppt%Tplgy(12)=INT(TmpTplgyArray(i-1,j+1,k  ))
        Ppt%Tplgy(13)=INT(TmpTplgyArray(i  ,j+1,k  ))
        Ppt%Tplgy(14)=INT(TmpTplgyArray(i+1,j+1,k  ))
        Ppt%Tplgy(15)=INT(TmpTplgyArray(i  ,j-1,k+1))
        Ppt%Tplgy(16)=INT(TmpTplgyArray(i-1,j  ,k+1))
        Ppt%Tplgy(17)=INT(TmpTplgyArray(i  ,j  ,k+1))
        Ppt%Tplgy(18)=INT(TmpTplgyArray(i+1,j  ,k+1))
        Ppt%Tplgy(19)=INT(TmpTplgyArray(i  ,j+1,k+1))
    ELSE
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "ERROR: Run_options_scheme command must be an integer between 1 and 2\n",ierr)
        STOP
    END IF
    CALL DMDAVecRestoreArrayReadF90(Gmtry%DstrMngr,Gmtry%Tplgy,TmpTplgyArray,ierr)

END SUBROUTINE GetLocalTopology

SUBROUTINE GeometryDestroy(Gmtry,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>

    PetscErrorCode,INTENT(INOUT)    :: ierr
    TYPE(Geometry),INTENT(IN)       :: Gmtry ! INOUT or IN?

    CALL DMDestroy(Gmtry%DstrMngr,ierr)
    CALL VecDestroy(Gmtry%Tplgy,ierr)

END SUBROUTINE GeometryDestroy

END MODULE ANISOFLOW_Geometry