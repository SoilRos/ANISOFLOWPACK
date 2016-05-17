MODULE ANISOFLOW_Properties

    USE ANISOFLOW_Interface, ONLY : GetVerbose
    USE ANISOFLOW_Types, ONLY : PropertyField,Geometry,Property

    IMPLICIT NONE

CONTAINS

SUBROUTINE GetProrperties(Gmtry,PptFld,ierr)
    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    TYPE(Geometry),INTENT(IN)           :: Gmtry
    TYPE(PropertyField),INTENT(OUT)     :: PptFld

    PetscLogStage                       :: Stage
    PetscBool                           :: Verbose

    CALL PetscLogStageRegister("GetProrperties", stage,ierr)
    CALL PetscLogStagePush(Stage,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetProrperties Stage] Inizialited\n",ierr)

    CALL GetConductivity(Gmtry,PptFld,ierr)

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetProrperties Stage] Finalized\n",ierr)
    CALL PetscLogStagePop(Stage,ierr)

END SUBROUTINE GetProrperties

SUBROUTINE GetConductivity(Gmtry,PptFld,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    TYPE(Geometry),INTENT(IN)           :: Gmtry
    TYPE(PropertyField),INTENT(INOUT)   :: PptFld

    TYPE(InputTypeVar)              :: InputType

    CALL GetInputType(InputType,ierr)

    IF (InputType%Tplgy.EQ.1) THEN
        CALL GetConductivity_1(Gmtry,PptFld,ierr)
    ELSE 
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "[ERROR] Conductivity InputType wrong\n",ierr)
        STOP
    END IF

END SUBROUTINE GetConductivity

SUBROUTINE GetConductivity_1(Gmtry,PptFld,ierr)

    USE ANISOFLOW_Interface
    USE ANISOFLOW_Types, ONLY : TargetFullTensor
    
    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscdmda.h90>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    TYPE(Geometry),INTENT(IN)           :: Gmtry
    TYPE(PropertyField),INTENT(INOUT)   :: PptFld

    PetscInt                            :: u,i,j,width(3),ValI,CvtLen
    PetscReal                           :: ValR
    PetscMPIInt                         :: process
    Vec                                 :: CvtTypeGlobal
    CHARACTER(LEN=200)                  :: InputDir,InputFileCvt,InputFileCvtByZones
    CHARACTER(LEN=200)                  :: Route
    CHARACTER(LEN=13)                   :: CvtKind
    PetscBool                           :: Verbose

    PARAMETER(u=01)

    CALL GetVerbose(Verbose,ierr)

    CALL GetInputDir(InputDir,ierr)
    CALL GetInputFileCvt(InputFileCvt,ierr)

    PptFld%Cvt%DefinedByZones=.TRUE.
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetProrperties Stage] Conductivity Field is stored as Zones by Block\n",ierr)

    CALL DMCreateGlobalVector(Gmtry%DataMngr,CvtTypeGlobal,ierr)
    CALL DMCreateLocalVector(Gmtry%DataMngr,PptFld%Cvt%CvtType,ierr)

    CALL MPI_Comm_rank(MPI_COMM_WORLD,process,ierr)
    IF (process.EQ.0) THEN
        CALL GetInputFileCvtByZones(InputFileCvtByZones,ierr)
        Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFileCvtByZones))
        OPEN(u,FILE=TRIM(Route),STATUS='OLD',ACTION='READ')
        READ(u, '((I10),(I10),(I10))')width(1),width(2),width(3)

        DO i=1,width(1)*width(2)*width(3)
            READ(u, '(I10)')ValI
            IF (ValI.LE.0) THEN
                CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,             &
                    & "[ERROR] Input_file_cvt_type entry only can contain"  &
                    & // " natural numbers\n",ierr)
                STOP
            END IF
            ValR=REAL(ValI)
            CALL VecSetValue(CvtTypeGlobal,i-1,ValR,INSERT_VALUES,ierr)
        END DO
        CLOSE(u)
    END IF

    CALL VecAssemblyBegin(CvtTypeGlobal,ierr)
    CALL VecAssemblyEnd(CvtTypeGlobal,ierr)

    CALL VecMax(CvtTypeGlobal,PETSC_NULL_INTEGER,ValR,ierr)
    CvtLen=INT(ValR)

    CALL DMGlobalToLocalBegin(Gmtry%DataMngr,CvtTypeGlobal,INSERT_VALUES,  &
        & PptFld%Cvt%CvtType,ierr)
    CALL DMGlobalToLocalEnd(Gmtry%DataMngr,CvtTypeGlobal,INSERT_VALUES,    &
        & PptFld%Cvt%CvtType,ierr)

    CALL VecDestroy(CvtTypeGlobal,ierr)

    ALLOCATE(PptFld%Cvt%CvtZone(CvtLen))
    
    IF (process.EQ.0) THEN
        Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFileCvt))
        OPEN(u,FILE=TRIM(Route),STATUS='OLD',ACTION='READ')

        DO i=1,CvtLen
            READ(u,*)
            READ(u,*)
            READ(u,*)
            READ(u,'(A13)')CvtKind
            IF (CvtKind.EQ."k anisotropic") THEN 
                READ(u,'((F8.0),(F8.0),(F8.0))')PptFld%Cvt%CvtZone(i)%xx, &
                    & PptFld%Cvt%CvtZone(i)%yy,PptFld%Cvt%CvtZone(i)%zz
                PptFld%Cvt%CvtZone(i)%xy=0.0
                PptFld%Cvt%CvtZone(i)%xz=0.0
                PptFld%Cvt%CvtZone(i)%yz=0.0
            ELSE IF (CvtKind.EQ."k isotropic  ") THEN
                READ(u,'(F15.0)')PptFld%Cvt%CvtZone(i)%xx
                PptFld%Cvt%CvtZone(i)%yy=PptFld%Cvt%CvtZone(i)%xx
                PptFld%Cvt%CvtZone(i)%zz=PptFld%Cvt%CvtZone(i)%xx
                PptFld%Cvt%CvtZone(i)%xy=0.0
                PptFld%Cvt%CvtZone(i)%xz=0.0
                PptFld%Cvt%CvtZone(i)%yz=0.0
            ELSE
                CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,             &
                    & "[ERROR] File of conductivity properties is invalid\n"&
                    & ,ierr)
                STOP            
            END IF
            DO j=1,24
                READ(u,*)
            END DO
        END DO
        CLOSE(u)
    END IF

    CALL MPI_Bcast(PptFld%Cvt%CvtZone(:)%xx,CvtLen,MPI_DOUBLE, 0,         &
        & PETSC_COMM_WORLD,ierr)
    CALL MPI_Bcast(PptFld%Cvt%CvtZone(:)%yy,CvtLen,MPI_DOUBLE, 0,         &
        & PETSC_COMM_WORLD,ierr)
    CALL MPI_Bcast(PptFld%Cvt%CvtZone(:)%zz,CvtLen,MPI_DOUBLE, 0,         &
        &PETSC_COMM_WORLD,ierr)
    CALL MPI_Bcast(PptFld%Cvt%CvtZone(:)%xy,CvtLen,MPI_DOUBLE, 0,         &
        &PETSC_COMM_WORLD,ierr)
    CALL MPI_Bcast(PptFld%Cvt%CvtZone(:)%xz,CvtLen,MPI_DOUBLE, 0,         &
        &PETSC_COMM_WORLD,ierr)
    CALL MPI_Bcast(PptFld%Cvt%CvtZone(:)%yz,CvtLen,MPI_DOUBLE, 0,         &
        &PETSC_COMM_WORLD,ierr)

    DO i=1,CvtLen
        CALL TargetFullTensor(PptFld%Cvt%CvtZone(i))
    END DO

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetProrperties Stage] Conductivity Field was satisfactorily created\n",ierr)

END SUBROUTINE GetConductivity_1

SUBROUTINE GetLocalProperty(Gmtry,PptFld,Ppt,i,j,k,ierr)

    USE ANISOFLOW_Geometry, ONLY : GetLocalTopology

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    TYPE(Geometry),INTENT(IN)           :: Gmtry
    TYPE(PropertyField),INTENT(IN)      :: PptFld
    TYPE(Property),INTENT(OUT)          :: Ppt
    PetscInt,INTENT(IN)                 :: i,j,k
    PetscErrorCode,INTENT(INOUT)        :: ierr

    Ppt%Pstn%i=i
    Ppt%Pstn%j=j
    Ppt%Pstn%k=k

    CALL GetLocalTopology(Gmtry,Ppt,ierr)
    CALL GetLocalConductivity(Gmtry,PptFld,Ppt,ierr)


END SUBROUTINE GetLocalProperty

SUBROUTINE GetLocalConductivity(Gmtry,PptFld,Ppt,ierr)

    USE ANISOFLOW_Interface
    USE ANISOFLOW_Types, ONLY : TargetFullTensor

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscdmda.h90>

    TYPE(Geometry),INTENT(IN)           :: Gmtry
    TYPE(PropertyField),INTENT(IN)      :: PptFld
    TYPE(Property),INTENT(INOUT)        :: Ppt
    PetscErrorCode,INTENT(INOUT)        :: ierr

    TYPE(InputTypeVar)                  :: InputType
    PetscReal,POINTER                   :: TmpCvtTypeZone(:,:,:)
    PetscInt                            :: ValI(2),i,j,k

    CALL GetInputType(InputType,ierr)

    i=Ppt%Pstn%i
    j=Ppt%Pstn%j
    k=Ppt%Pstn%k

    Ppt%CvtBx%xx=0.0
    Ppt%CvtBx%xy=0.0
    Ppt%CvtBx%xz=0.0
    Ppt%CvtBx%yy=0.0
    Ppt%CvtBx%yz=0.0
    Ppt%CvtBx%zz=0.0

    Ppt%CvtFx%xx=0.0
    Ppt%CvtFx%xy=0.0
    Ppt%CvtFx%xz=0.0
    Ppt%CvtFx%yy=0.0
    Ppt%CvtFx%yz=0.0
    Ppt%CvtFx%zz=0.0

    Ppt%CvtBy%xx=0.0
    Ppt%CvtBy%xy=0.0
    Ppt%CvtBy%xz=0.0
    Ppt%CvtBy%yy=0.0
    Ppt%CvtBy%yz=0.0
    Ppt%CvtBy%zz=0.0

    Ppt%CvtFy%xx=0.0
    Ppt%CvtFy%xy=0.0
    Ppt%CvtFy%xz=0.0
    Ppt%CvtFy%yy=0.0
    Ppt%CvtFy%yz=0.0
    Ppt%CvtFy%zz=0.0

    Ppt%CvtBz%xx=0.0
    Ppt%CvtBz%xy=0.0
    Ppt%CvtBz%xz=0.0
    Ppt%CvtBz%yy=0.0
    Ppt%CvtBz%yz=0.0
    Ppt%CvtBz%zz=0.0

    Ppt%CvtFz%xx=0.0
    Ppt%CvtFz%xy=0.0
    Ppt%CvtFz%xz=0.0
    Ppt%CvtFz%yy=0.0
    Ppt%CvtFz%yz=0.0
    Ppt%CvtFz%zz=0.0

    IF (SIZE(Ppt%StnclTplgy).EQ.7) THEN          ! Stencil type star
        ValI(1)=4
    ELSEIF (SIZE(Ppt%StnclTplgy).EQ.19) THEN     ! Stencil type box
        ValI(1)=10
    ELSE 
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,              &
            & "[ERROR] Property topology wasn't well defined\n"      &
            & ,ierr)
        STOP   
    END IF

    IF (.NOT.PptFld%Cvt%DefinedByInterface) THEN
        IF ((Ppt%StnclTplgy(ValI(1)).EQ.1).OR.(Ppt%StnclTplgy(ValI(1)).EQ.3).OR.(Ppt%StnclTplgy(ValI(1)).EQ.4).OR.(Ppt%StnclTplgy(ValI(1)).EQ.5)) THEN ! Only get properties for active blocks
            CALL DMDAVecGetArrayreadF90(Gmtry%DataMngr,PptFld%Cvt%CvtType,     &
                & TmpCvtTypeZone,ierr)

            ValI(1)=INT(TmpCvtTypeZone(i,j,k))

            ! Backward on x
            ValI(2)=INT(TmpCvtTypeZone(i-1,j,k))
            Ppt%CvtBx%xx=Armonic(PptFld%Cvt%CvtZone(ValI(1))%xx,PptFld%Cvt%CvtZone(ValI(2))%xx)
            Ppt%CvtBx%yy=Armonic(PptFld%Cvt%CvtZone(ValI(1))%yy,PptFld%Cvt%CvtZone(ValI(2))%yy)
            Ppt%CvtBx%zz=Armonic(PptFld%Cvt%CvtZone(ValI(1))%zz,PptFld%Cvt%CvtZone(ValI(2))%zz)
            Ppt%CvtBx%xy=Armonic(PptFld%Cvt%CvtZone(ValI(1))%xy,PptFld%Cvt%CvtZone(ValI(2))%xy)
            Ppt%CvtBx%xz=Armonic(PptFld%Cvt%CvtZone(ValI(1))%xz,PptFld%Cvt%CvtZone(ValI(2))%xz)
            Ppt%CvtBx%yz=Armonic(PptFld%Cvt%CvtZone(ValI(1))%yz,PptFld%Cvt%CvtZone(ValI(2))%yz)
            CALL TargetFullTensor(Ppt%CvtBx)

            ! Forward on x
            ValI(2)=INT(TmpCvtTypeZone(i+1,j,k))
            Ppt%CvtFx%xx=Armonic(PptFld%Cvt%CvtZone(ValI(1))%xx,PptFld%Cvt%CvtZone(ValI(2))%xx)
            Ppt%CvtFx%yy=Armonic(PptFld%Cvt%CvtZone(ValI(1))%yy,PptFld%Cvt%CvtZone(ValI(2))%yy)
            Ppt%CvtFx%zz=Armonic(PptFld%Cvt%CvtZone(ValI(1))%zz,PptFld%Cvt%CvtZone(ValI(2))%zz)
            Ppt%CvtFx%xy=Armonic(PptFld%Cvt%CvtZone(ValI(1))%xy,PptFld%Cvt%CvtZone(ValI(2))%xy)
            Ppt%CvtFx%xz=Armonic(PptFld%Cvt%CvtZone(ValI(1))%xz,PptFld%Cvt%CvtZone(ValI(2))%xz)
            Ppt%CvtFx%yz=Armonic(PptFld%Cvt%CvtZone(ValI(1))%yz,PptFld%Cvt%CvtZone(ValI(2))%yz)
            CALL TargetFullTensor(Ppt%CvtFx)

            ! Backward on y
            ValI(2)=INT(TmpCvtTypeZone(i,j-1,k))
            Ppt%CvtBy%xx=Armonic(PptFld%Cvt%CvtZone(ValI(1))%xx,PptFld%Cvt%CvtZone(ValI(2))%xx)
            Ppt%CvtBy%yy=Armonic(PptFld%Cvt%CvtZone(ValI(1))%yy,PptFld%Cvt%CvtZone(ValI(2))%yy)
            Ppt%CvtBy%zz=Armonic(PptFld%Cvt%CvtZone(ValI(1))%zz,PptFld%Cvt%CvtZone(ValI(2))%zz)
            Ppt%CvtBy%xy=Armonic(PptFld%Cvt%CvtZone(ValI(1))%xy,PptFld%Cvt%CvtZone(ValI(2))%xy)
            Ppt%CvtBy%xz=Armonic(PptFld%Cvt%CvtZone(ValI(1))%xz,PptFld%Cvt%CvtZone(ValI(2))%xz)
            Ppt%CvtBy%yz=Armonic(PptFld%Cvt%CvtZone(ValI(1))%yz,PptFld%Cvt%CvtZone(ValI(2))%yz)
            CALL TargetFullTensor(Ppt%CvtBy)

            ! Forward on y
            ValI(2)=INT(TmpCvtTypeZone(i,j+1,k))
            Ppt%CvtFy%xx=Armonic(PptFld%Cvt%CvtZone(ValI(1))%xx,PptFld%Cvt%CvtZone(ValI(2))%xx)
            Ppt%CvtFy%yy=Armonic(PptFld%Cvt%CvtZone(ValI(1))%yy,PptFld%Cvt%CvtZone(ValI(2))%yy)
            Ppt%CvtFy%zz=Armonic(PptFld%Cvt%CvtZone(ValI(1))%zz,PptFld%Cvt%CvtZone(ValI(2))%zz)
            Ppt%CvtFy%xy=Armonic(PptFld%Cvt%CvtZone(ValI(1))%xy,PptFld%Cvt%CvtZone(ValI(2))%xy)
            Ppt%CvtFy%xz=Armonic(PptFld%Cvt%CvtZone(ValI(1))%xz,PptFld%Cvt%CvtZone(ValI(2))%xz)
            Ppt%CvtFy%yz=Armonic(PptFld%Cvt%CvtZone(ValI(1))%yz,PptFld%Cvt%CvtZone(ValI(2))%yz)
            CALL TargetFullTensor(Ppt%CvtFy)

            ! Backward on z
            ValI(2)=INT(TmpCvtTypeZone(i,j,k-1))
            Ppt%CvtBz%xx=Armonic(PptFld%Cvt%CvtZone(ValI(1))%xx,PptFld%Cvt%CvtZone(ValI(2))%xx)
            Ppt%CvtBz%yy=Armonic(PptFld%Cvt%CvtZone(ValI(1))%yy,PptFld%Cvt%CvtZone(ValI(2))%yy)
            Ppt%CvtBz%zz=Armonic(PptFld%Cvt%CvtZone(ValI(1))%zz,PptFld%Cvt%CvtZone(ValI(2))%zz)
            Ppt%CvtBz%xy=Armonic(PptFld%Cvt%CvtZone(ValI(1))%xy,PptFld%Cvt%CvtZone(ValI(2))%xy)
            Ppt%CvtBz%xz=Armonic(PptFld%Cvt%CvtZone(ValI(1))%xz,PptFld%Cvt%CvtZone(ValI(2))%xz)
            Ppt%CvtBz%yz=Armonic(PptFld%Cvt%CvtZone(ValI(1))%yz,PptFld%Cvt%CvtZone(ValI(2))%yz)
            CALL TargetFullTensor(Ppt%CvtBz)

            ! Forward on y
            ValI(2)=INT(TmpCvtTypeZone(i,j,k+1))
            Ppt%CvtFz%xx=Armonic(PptFld%Cvt%CvtZone(ValI(1))%xx,PptFld%Cvt%CvtZone(ValI(2))%xx)
            Ppt%CvtFz%yy=Armonic(PptFld%Cvt%CvtZone(ValI(1))%yy,PptFld%Cvt%CvtZone(ValI(2))%yy)
            Ppt%CvtFz%zz=Armonic(PptFld%Cvt%CvtZone(ValI(1))%zz,PptFld%Cvt%CvtZone(ValI(2))%zz)
            Ppt%CvtFz%xy=Armonic(PptFld%Cvt%CvtZone(ValI(1))%xy,PptFld%Cvt%CvtZone(ValI(2))%xy)
            Ppt%CvtFz%xz=Armonic(PptFld%Cvt%CvtZone(ValI(1))%xz,PptFld%Cvt%CvtZone(ValI(2))%xz)
            Ppt%CvtFz%yz=Armonic(PptFld%Cvt%CvtZone(ValI(1))%yz,PptFld%Cvt%CvtZone(ValI(2))%yz)
            CALL TargetFullTensor(Ppt%CvtFz)

            CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%CvtType, &
                & TmpCvtTypeZone,ierr)

        END IF
    ELSE
          ! TO DO: It depends on how we're going to save the Interface conductivity
    END IF

END SUBROUTINE GetLocalConductivity

PetscReal FUNCTION Armonic(ValR1,ValR2)
    IMPLICIT NONE
#include <petsc/finclude/petsc.h>
    PetscReal, INTENT(IN) :: ValR1,ValR2
    IF ((ValR1==0).OR.(ValR2==0)) THEN
        Armonic=0.0
    ELSE
        Armonic = 2.0/(1.0/ValR1+1.0/ValR2)
    END IF
END FUNCTION Armonic

SUBROUTINE DestroyProperties(PptFld,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    TYPE(PropertyField),INTENT(INOUT)   :: PptFld

    TYPE(InputTypeVar)                  :: InputType
    PetscLogStage                       :: Stage
    PetscBool                           :: Verbose

    CALL PetscLogStageRegister("DestroyProperties", stage,ierr)
    CALL PetscLogStagePush(Stage,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[DestroyProperties Stage] Inizialited\n",ierr)

    CALL GetInputType(InputType,ierr)

    IF (PptFld%Cvt%DefinedByZones) THEN
        IF (ALLOCATED(PptFld%Cvt%CvtZone)) DEALLOCATE(PptFld%Cvt%CvtZone)
        CALL VecDestroy(PptFld%Cvt%CvtType,ierr)
    ELSE
        CALL VecDestroy(PptFld%Cvt%xxVec,ierr)
        CALL VecDestroy(PptFld%Cvt%yyVec,ierr)
        CALL VecDestroy(PptFld%Cvt%zzVec,ierr)
    END IF

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[DestroyProperties Stage] Finalized\n",ierr)
    CALL PetscLogStagePop(Stage,ierr)

END SUBROUTINE DestroyProperties

END MODULE ANISOFLOW_Properties