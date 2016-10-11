MODULE ANISOFLOW_Properties

    IMPLICIT NONE

CONTAINS

SUBROUTINE GetProrperties(Gmtry,PptFld,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,PropertiesField,RunOptionsVar
    USE ANISOFLOW_Interface, ONLY : GetVerbose,GetRunOptions

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscdmda.h90>


    PetscErrorCode,INTENT(INOUT)        :: ierr
    TYPE(Geometry),INTENT(IN)           :: Gmtry
    TYPE(PropertiesField),INTENT(OUT)   :: PptFld

    CHARACTER(LEN=200)                  :: EventName,ClassName
    PetscBool                           :: Verbose
    PetscLogEvent                       :: Event
    PetscClassId                        :: ClassID
    PetscLogDouble                      :: EventFlops=0.D0
    TYPE(RunOptionsVar)                 :: RunOptions

    ClassName="Property"
    CALL PetscClassIdRegister(ClassName,ClassID,ierr)
    EventName="GetProrperties"
    CALL PetscLogEventRegister(EventName,ClassID,Event,ierr)
    CALL PetscLogEventBegin(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Inizialited\n",ierr)

    CALL GetPptZoneID(Gmtry,PptFld%ZoneID,PptFld%DefinedByPptZones,ierr)
    CALL GetConductivity(Gmtry,PptFld,PptFld%Cvt,ierr)

    CALL GetRunOptions(RunOptions,ierr)
    IF (RunOptions%Time) THEN ! Transitory
        CALL GetStorage(Gmtry,PptFld,PptFld%Sto,ierr)
    END IF

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Finalized\n",ierr)
    
    CALL PetscLogFlops(EventFlops,ierr)
    CALL PetscLogEventEnd(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)


END SUBROUTINE GetProrperties

SUBROUTINE GetPptZoneID(Gmtry,PptZoneID_Local,DefinedByPptZones,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry
    USE ANISOFLOW_Interface, ONLY : GetVerbose,GetInputDir,GetInputFilePptByZone,GetHomogeneusPptFlg
    USE ANISOFLOW_View, ONLY : ViewProperty
    USE ANISOFLOW_Geometry, ONLY : VecApplicationToPetsc

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscdmda.h90>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    TYPE(Geometry),INTENT(IN)           :: Gmtry
    Vec,INTENT(OUT)                     :: PptZoneID_Local
    PetscBool,INTENT(OUT)               :: DefinedByPptZones

    PetscInt                            :: widthG(3),u,ValI,i
    PetscMPIInt                         :: process
    PetscReal                           :: ValR,one=1.D0
    CHARACTER(LEN=200)                  :: InputFilePptByZone,Route,ViewName,InputDir,EventName
    PetscBool                           :: InputFilePptByZoneFlg,HomogeneusPptFlg,Verbose
    Vec                                 :: PptZoneID_Global

    PARAMETER(u=01)

    CALL GetVerbose(Verbose,ierr)

    CALL GetInputDir(InputDir,ierr)
    CALL GetInputFilePptByZone(InputFilePptByZone,InputFilePptByZoneFlg,ierr)
    CALL GetHomogeneusPptFlg(HomogeneusPptFlg,ierr)

    IF (InputFilePptByZoneFlg.AND.(.NOT.HomogeneusPptFlg)) THEN

        CALL DMCreateGlobalVector(Gmtry%DataMngr,PptZoneID_Global,ierr)
        CALL DMCreateLocalVector(Gmtry%DataMngr,PptZoneID_Local,ierr)

        CALL DMDAGetInfo(Gmtry%DataMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
            & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
            & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
            & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
            & PETSC_NULL_INTEGER,ierr)

        CALL MPI_Comm_rank(MPI_COMM_WORLD,process,ierr)


            IF (process.EQ.0) THEN

                Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFilePptByZone))
                OPEN(u,FILE=TRIM(Route),STATUS='OLD',ACTION='READ')

                DO i=1,widthG(1)*widthG(2)*widthG(3)
                    READ(u,*)ValI
                    IF (ValI.LE.0) THEN
                        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,             &
                            & "[ERROR] -Input_file_ppt_by_zones entry only can contain"  &
                            & // " natural numbers\n",ierr)
                        STOP
                    END IF
                    ValR=REAL(ValI)
                    CALL VecSetValue(PptZoneID_Global,i-1,ValR,INSERT_VALUES,ierr)
                END DO
                CLOSE(u)

            END IF

        CALL VecAssemblyBegin(PptZoneID_Global,ierr)
        CALL VecAssemblyEnd(PptZoneID_Global,ierr)

        CALL VecApplicationToPetsc(Gmtry%DataMngr,PptZoneID_Global,ierr)

        ViewName="PptZoneID"
        EventName="GetPptZoneID"
        CALL ViewProperty(PptZoneID_Global,ViewName,EventName,ierr)

        CALL DMGlobalToLocalBegin(Gmtry%DataMngr,PptZoneID_Global,INSERT_VALUES,PptZoneID_Local,ierr)
        CALL DMGlobalToLocalEnd(Gmtry%DataMngr,PptZoneID_Global,INSERT_VALUES,PptZoneID_Local,ierr)
        CALL VecDestroy(PptZoneID_Global,ierr)
        DefinedByPptZones=.TRUE.

    ELSEIF (HomogeneusPptFlg) THEN

        CALL DMCreateGlobalVector(Gmtry%DataMngr,PptZoneID_Global,ierr)
        CALL DMCreateLocalVector(Gmtry%DataMngr,PptZoneID_Local,ierr)

        CALL DMDAGetInfo(Gmtry%DataMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
            & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
            & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
            & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
            & PETSC_NULL_INTEGER,ierr)

        CALL MPI_Comm_rank(MPI_COMM_WORLD,process,ierr)

        CALL VecSet(PptZoneID_Global,one,ierr)

        CALL VecAssemblyBegin(PptZoneID_Global,ierr)
        CALL VecAssemblyEnd(PptZoneID_Global,ierr)

        CALL VecApplicationToPetsc(Gmtry%DataMngr,PptZoneID_Global,ierr)

        ViewName="PptZoneID"
        EventName="GetPptZoneID"
        CALL ViewProperty(PptZoneID_Global,ViewName,EventName,ierr)

        CALL DMGlobalToLocalBegin(Gmtry%DataMngr,PptZoneID_Global,INSERT_VALUES,PptZoneID_Local,ierr)
        CALL DMGlobalToLocalEnd(Gmtry%DataMngr,PptZoneID_Global,INSERT_VALUES,PptZoneID_Local,ierr)
        CALL VecDestroy(PptZoneID_Global,ierr)
        DefinedByPptZones=.TRUE.
        IF (InputFilePptByZoneFlg) THEN 
            IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetProrperties Event] WARNING: -Input_file_ppt_by_zones and -Homogeneuos_ppt used at the same time; properties was set homogeneuos in the whole domain..\n",ierr)
        END IF
    END IF

END SUBROUTINE GetPptZoneID

SUBROUTINE GetConductivity(Gmtry,PptFld,Cvt,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,PropertiesField,ConductivityField,InputTypeVar
    USE ANISOFLOW_Interface, ONLY : GetVerbose,GetInputType

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    TYPE(Geometry),INTENT(IN)           :: Gmtry
    TYPE(PropertiesField),INTENT(IN)    :: PptFld
    TYPE(ConductivityField),INTENT(OUT) :: Cvt

    TYPE(InputTypeVar)                  :: InputType
    PetscBool                           :: Verbose
    CHARACTER(LEN=200)                  :: EventName,ClassName
    PetscLogEvent                       :: Event
    PetscClassId                        :: ClassID
    PetscLogDouble                      :: EventFlops=0.d0

    ClassName="Property"
    CALL PetscClassIdRegister(ClassName,ClassID,ierr)
    EventName="GetConductivity"
    CALL PetscLogEventRegister(EventName,ClassID,Event,ierr)
    CALL PetscLogEventBegin(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Inizialited\n",ierr)
    
    CALL GetInputType(InputType,ierr)

    IF (InputType%Cvt.EQ.1) THEN
        CALL GetConductivity_1(Gmtry,PptFld,Cvt,ierr)
    ELSE IF (InputType%Cvt.EQ.2) THEN
        CALL GetConductivity_2(Gmtry,PptFld,Cvt,ierr)
    ELSE IF (InputType%Cvt.EQ.3) THEN
        CALL GetConductivity_3(Gmtry,PptFld,Cvt,ierr)
    ELSE IF (InputType%Cvt.EQ.4) THEN
        CALL GetConductivity_4(Gmtry,PptFld,Cvt,ierr)
    ELSE
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "[ERROR] Conductivity InputType wrong\n",ierr)
        STOP
    END IF

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Finalized\n",ierr)
    
    CALL PetscLogFlops(EventFlops,ierr)
    CALL PetscLogEventEnd(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

END SUBROUTINE GetConductivity

SUBROUTINE GetCvtZoneID(Gmtry,PptFld,CvtZoneID_Local,DefinedBy,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,PropertiesField,TargPetscVec
    USE ANISOFLOW_Interface, ONLY : GetVerbose,GetInputDir,GetInputFilePptByZone,GetInputFileCvtByZone,GetHomogeneusPptFlg,GetHomogeneusCvtFlg
    USE ANISOFLOW_View, ONLY : ViewConductivity
    USE ANISOFLOW_Geometry, ONLY : VecApplicationToPetsc

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscdmda.h90>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    TYPE(Geometry),INTENT(IN)           :: Gmtry
    TYPE(PropertiesField),INTENT(IN)    :: PptFld
    Vec,INTENT(OUT)                     :: CvtZoneID_Local
    PetscInt,INTENT(OUT)                :: DefinedBy

    PetscInt                            :: widthG(3),u,ValI,i
    PetscMPIInt                         :: process
    PetscReal                           :: ValR,one=1.D0
    CHARACTER(LEN=200)                  :: InputFilePptByZone,InputFileCvtByZone,Route,ViewName,InputDir,EventName
    PetscBool                           :: InputFilePptByZoneFlg,InputFileCvtByZoneFlg,HomogeneusPptFlg,HomogeneusCvtFlg,Verbose
    Vec                                 :: CvtZoneID_Global
    Vec,POINTER                         :: ZoneID_tmp

    PARAMETER(u=01)

    CALL GetVerbose(Verbose,ierr)

    CALL GetInputDir(InputDir,ierr)
    CALL GetInputFilePptByZone(InputFilePptByZone,InputFilePptByZoneFlg,ierr)
    CALL GetInputFileCvtByZone(InputFileCvtByZone,InputFileCvtByZoneFlg,ierr)

    CALL GetHomogeneusCvtFlg(HomogeneusCvtFlg,ierr)
    CALL GetHomogeneusPptFlg(HomogeneusPptFlg,ierr)

    IF (InputFileCvtByZoneFlg.AND.(.NOT.HomogeneusCvtFlg)) THEN

        CALL DMCreateGlobalVector(Gmtry%DataMngr,CvtZoneID_Global,ierr)
        CALL DMCreateLocalVector(Gmtry%DataMngr,CvtZoneID_Local,ierr)

        CALL DMDAGetInfo(Gmtry%DataMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
            & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
            & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
            & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
            & PETSC_NULL_INTEGER,ierr)

        CALL MPI_Comm_rank(MPI_COMM_WORLD,process,ierr)
            IF (process.EQ.0) THEN
                Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFileCvtByZone))
                OPEN(u,FILE=TRIM(Route),STATUS='OLD',ACTION='READ')
                DO i=1,widthG(1)*widthG(2)*widthG(3)
                    READ(u,*)ValI
                    IF (ValI.LE.0) THEN
                        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,             &
                            & "[ERROR] Input_file_cvt_by_zones entry only can contain"  &
                            & // " natural numbers\n",ierr)
                        STOP
                    END IF
                    ValR=REAL(ValI)
                    CALL VecSetValue(CvtZoneID_Global,i-1,ValR,INSERT_VALUES,ierr)
                END DO
                CLOSE(u)
            END IF

        CALL VecAssemblyBegin(CvtZoneID_Global,ierr)
        CALL VecAssemblyEnd(CvtZoneID_Global,ierr)

        CALL VecApplicationToPetsc(Gmtry%DataMngr,CvtZoneID_Global,ierr)

        ViewName="CvtZoneID"
        EventName="GetCvtZoneID"
        CALL ViewConductivity(CvtZoneID_Global,ViewName,EventName,ierr)

        CALL DMGlobalToLocalBegin(Gmtry%DataMngr,CvtZoneID_Global,INSERT_VALUES,CvtZoneID_Local,ierr)
        CALL DMGlobalToLocalEnd(Gmtry%DataMngr,CvtZoneID_Global,INSERT_VALUES,CvtZoneID_Local,ierr)
        CALL VecDestroy(CvtZoneID_Global,ierr)

        DefinedBy=1
    ELSEIF (HomogeneusCvtFlg) THEN

        CALL DMCreateGlobalVector(Gmtry%DataMngr,CvtZoneID_Global,ierr)
        CALL DMCreateLocalVector(Gmtry%DataMngr,CvtZoneID_Local,ierr)

        CALL DMDAGetInfo(Gmtry%DataMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
            & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
            & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
            & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
            & PETSC_NULL_INTEGER,ierr)

        CALL MPI_Comm_rank(MPI_COMM_WORLD,process,ierr)

        CALL VecSet(CvtZoneID_Global,one,ierr)

        CALL VecAssemblyBegin(CvtZoneID_Global,ierr)
        CALL VecAssemblyEnd(CvtZoneID_Global,ierr)

        CALL VecApplicationToPetsc(Gmtry%DataMngr,CvtZoneID_Global,ierr)

        ViewName="CvtZoneID"
        EventName="GetCvtZoneID"
        CALL ViewConductivity(CvtZoneID_Global,ViewName,EventName,ierr)

        CALL DMGlobalToLocalBegin(Gmtry%DataMngr,CvtZoneID_Global,INSERT_VALUES,CvtZoneID_Local,ierr)
        CALL DMGlobalToLocalEnd(Gmtry%DataMngr,CvtZoneID_Global,INSERT_VALUES,CvtZoneID_Local,ierr)
        CALL VecDestroy(CvtZoneID_Global,ierr)

        DefinedBy=1
        IF (InputFileCvtByZoneFlg) THEN 
            IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetProrperties Event] WARNING: -Input_file_cvt_by_zones and -Homogeneuos_cvt used at the same time; conductivities was set homogeneuos in the whole domain.\n",ierr)
        END IF
    ELSEIF (InputFilePptByZoneFlg.OR.HomogeneusPptFlg) THEN
        ZoneID_tmp => TargPetscVec(PptFld%ZoneID)
        CvtZoneID_Local = ZoneID_tmp
        DefinedBy=2
    ELSE
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,             &
            & "[ERROR] -Input_type_cvt 1 must be used with -Input_file_(cvt,ppt)_by_zones,"  &
            & // " which contains the conductivities zonification, or with -Homogeneuos_(cvt,ppt) 1, which sets the first property of -Input_file_(cvt,ppt) to the whole domain.\n",ierr)
        STOP
    END IF

END SUBROUTINE GetCvtZoneID

SUBROUTINE GetConductivity_1(Gmtry,PptFld,Cvt,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,ConductivityField,PropertiesField,TargetFullTensor
    USE ANISOFLOW_Interface, ONLY : GetVerbose,GetInputDir,GetInputFileCvt
    USE ANISOFLOW_Geometry, ONLY : VecApplicationToPetsc
    
    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscdmda.h90>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    TYPE(Geometry),INTENT(IN)           :: Gmtry
    TYPE(PropertiesField),INTENT(IN)    :: PptFld
    TYPE(ConductivityField),INTENT(OUT) :: Cvt

    PetscInt                            :: u,i,j,CvtLen
    PetscReal                           :: ValR
    PetscMPIInt                         :: process
    CHARACTER(LEN=200)                  :: InputDir,InputFileCvt
    CHARACTER(LEN=200)                  :: Route
    CHARACTER(LEN=200)                  :: CvtKind
    PetscBool                           :: Verbose

    PARAMETER(u=01)

    CALL GetVerbose(Verbose,ierr)

    CALL GetCvtZoneID(Gmtry,PptFld,Cvt%ZoneID,Cvt%DefinedBy,ierr)

    CALL GetInputDir(InputDir,ierr)
    CALL GetInputFileCvt(InputFileCvt,ierr)

    CALL MPI_Comm_rank(MPI_COMM_WORLD,process,ierr)

!     IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetConductivity Event] Conductivity Field is stored as Zones by Block\n",ierr)

    CALL VecMax(Cvt%ZoneID,PETSC_NULL_INTEGER,ValR,ierr)
    CvtLen=INT(ValR)

    ALLOCATE(Cvt%Zone(CvtLen))

    IF (process.EQ.0) THEN
        Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFileCvt))
        OPEN(u,FILE=TRIM(Route),STATUS='OLD',ACTION='READ')

        DO i=1,CvtLen
            READ(u,*)
            READ(u,*)
            READ(u,*)
            READ(u,'(A)')CvtKind

            IF (ADJUSTL(TRIM(CvtKind)).EQ."k isotropic  ") THEN
                READ(u,*)Cvt%Zone(i)%xx
                Cvt%Zone(i)%yy=Cvt%Zone(i)%xx
                Cvt%Zone(i)%zz=Cvt%Zone(i)%xx
                Cvt%Zone(i)%xy=0.D0
                Cvt%Zone(i)%xz=0.D0
                Cvt%Zone(i)%yz=0.D0
            ELSEIF (ADJUSTL(TRIM(CvtKind)).EQ."k anisotropic") THEN 
                READ(u,*)Cvt%Zone(i)%xx, &
                    & Cvt%Zone(i)%yy,Cvt%Zone(i)%zz
                Cvt%Zone(i)%xy=0.D0
                Cvt%Zone(i)%xz=0.D0
                Cvt%Zone(i)%yz=0.D0
            ELSEIF (ADJUSTL(TRIM(CvtKind)).EQ."k non-orthogonal anisotropic") THEN 
                READ(u,*)Cvt%Zone(i)%xx,Cvt%Zone(i)%yy,Cvt%Zone(i)%zz,&
                        &Cvt%Zone(i)%xy,Cvt%Zone(i)%xz,Cvt%Zone(i)%yz
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

    CALL MPI_Bcast(Cvt%Zone(:)%xx,CvtLen,MPI_DOUBLE,0,PETSC_COMM_WORLD,ierr)
    CALL MPI_Bcast(Cvt%Zone(:)%yy,CvtLen,MPI_DOUBLE,0,PETSC_COMM_WORLD,ierr)
    CALL MPI_Bcast(Cvt%Zone(:)%zz,CvtLen,MPI_DOUBLE,0,PETSC_COMM_WORLD,ierr)
    CALL MPI_Bcast(Cvt%Zone(:)%xy,CvtLen,MPI_DOUBLE,0,PETSC_COMM_WORLD,ierr)
    CALL MPI_Bcast(Cvt%Zone(:)%xz,CvtLen,MPI_DOUBLE,0,PETSC_COMM_WORLD,ierr)
    CALL MPI_Bcast(Cvt%Zone(:)%yz,CvtLen,MPI_DOUBLE,0,PETSC_COMM_WORLD,ierr)

    DO i=1,CvtLen
        CALL TargetFullTensor(Cvt%Zone(i))
    END DO

!     IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetConductivity Event] Conductivity Field was satisfactorily created\n",ierr)

END SUBROUTINE GetConductivity_1

SUBROUTINE GetConductivity_2(Gmtry,PptFld,Cvt,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,PropertiesField,ConductivityField,TargPetscVec
    USE ANISOFLOW_Interface, ONLY : GetVerbose,GetInputDir,GetInputFileCvt
    USE ANISOFLOW_View, ONLY : ViewConductivity
    
    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscdmda.h90>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    TYPE(Geometry),INTENT(IN)           :: Gmtry
    TYPE(PropertiesField),INTENT(IN)    :: PptFld
    TYPE(ConductivityField),INTENT(OUT) :: Cvt

    PetscInt                            :: u,i,widthG(3)
    PetscReal                           :: ValR
    PetscMPIInt                         :: process
    CHARACTER(LEN=200)                  :: InputDir,InputFileCvt
    CHARACTER(LEN=200)                  :: Route,ViewName,EventName
    PetscBool                           :: Verbose
    Vec,POINTER                         :: Cvt_yy_pointer,Cvt_zz_pointer,Cvt_xz_pointer,Cvt_yz_pointer
    Vec                                 :: Cvt_xx_Global

    PARAMETER(u=01)

    CALL GetVerbose(Verbose,ierr)

    CALL GetInputDir(InputDir,ierr)
    CALL GetInputFileCvt(InputFileCvt,ierr)

    Cvt%DefinedBy=3
!     IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetProrperties Event] Conductivity Field is stored by Block\n",ierr)

    CALL DMCreateLocalVector(Gmtry%DataMngr,Cvt%xx,ierr)
    CALL DMCreateGlobalVector(Gmtry%DataMngr,Cvt_xx_Global,ierr)

    CALL MPI_Comm_rank(MPI_COMM_WORLD,process,ierr)

    IF (process.EQ.0) THEN
        Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFileCvt))
        OPEN(u,FILE=TRIM(Route),STATUS='OLD',ACTION='READ')

        ! It gets the global size from the geometry data manager.
        CALL DMDAGetInfo(Gmtry%DataMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
        & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,ierr)

        DO i=1,widthG(1)*widthG(2)*widthG(3)
            READ(u,*)ValR
            CALL VecSetValue(Cvt_xx_Global,i-1,ValR,INSERT_VALUES,ierr)
        END DO
        CLOSE(u)
    END IF

    CALL VecAssemblyBegin(Cvt_xx_Global,ierr)
    CALL VecAssemblyEnd(Cvt_xx_Global,ierr)

    ViewName="Cvt"
    EventName="GetConductivity"
    CALL ViewConductivity(Cvt%xx,ViewName,EventName,ierr)

    CALL DMGlobalToLocalBegin(Gmtry%DataMngr,Cvt_xx_Global,INSERT_VALUES,Cvt%xx,ierr)
    CALL DMGlobalToLocalEnd(Gmtry%DataMngr,Cvt_xx_Global,INSERT_VALUES,Cvt%xx,ierr)
    CALL VecDestroy(Cvt_xx_Global,ierr)

    Cvt_yy_pointer => TargPetscVec(Cvt%xx)
    Cvt%yy = Cvt_yy_pointer
    Cvt_zz_pointer => TargPetscVec(Cvt%xx)
    Cvt%zz = Cvt_zz_pointer

    ! Saving a vector with zeros (If we put a conditionals on the stencil it can be avoided, but would be unclear)
    CALL DMCreateLocalVector(Gmtry%DataMngr,Cvt%xy,ierr)
    CALL VecZeroEntries(Cvt%xy,ierr)
    CALL VecAssemblyBegin(Cvt%xy,ierr)
    CALL VecAssemblyEnd(Cvt%xy,ierr)

    Cvt_xz_pointer => TargPetscVec(Cvt%xy)
    Cvt%xz = Cvt_xz_pointer
    Cvt_yz_pointer => TargPetscVec(Cvt%xy)
    Cvt%yz = Cvt_yz_pointer

!     IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetProrperties Stage] Conductivity Field was satisfactorily created\n",ierr)

END SUBROUTINE GetConductivity_2

SUBROUTINE GetConductivity_3(Gmtry,PptFld,Cvt,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,PropertiesField,ConductivityField,TargPetscVec
    USE ANISOFLOW_Interface, ONLY : GetVerbose,GetInputDir,GetInputFileCvt
    USE ANISOFLOW_View, ONLY : ViewConductivity
    
    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscdmda.h90>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    TYPE(Geometry),INTENT(IN)           :: Gmtry
    TYPE(PropertiesField),INTENT(IN)    :: PptFld
    TYPE(ConductivityField),INTENT(OUT) :: Cvt

    PetscInt                            :: u,i,widthG(3)
    PetscReal                           :: ValR(3)
    PetscMPIInt                         :: process
    CHARACTER(LEN=200)                  :: InputDir,InputFileCvt
    CHARACTER(LEN=200)                  :: Route,ViewName,EventName
    PetscBool                           :: Verbose
    Vec,POINTER                         :: Cvt_xz_pointer,Cvt_yz_pointer
    Vec                                 :: Cvt_xx_Global,Cvt_yy_Global,Cvt_zz_Global

    PARAMETER(u=01)

    CALL GetVerbose(Verbose,ierr)

    CALL GetInputDir(InputDir,ierr)
    CALL GetInputFileCvt(InputFileCvt,ierr)

    Cvt%DefinedBy=3
!     IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetProrperties Event] Conductivity Field is stored by Block\n",ierr)

    CALL DMCreateLocalVector(Gmtry%DataMngr,Cvt%xx,ierr)
    CALL DMCreateLocalVector(Gmtry%DataMngr,Cvt%yy,ierr)
    CALL DMCreateLocalVector(Gmtry%DataMngr,Cvt%zz,ierr)

    CALL DMCreateGlobalVector(Gmtry%DataMngr,Cvt_xx_Global,ierr)
    CALL DMCreateGlobalVector(Gmtry%DataMngr,Cvt_yy_Global,ierr)
    CALL DMCreateGlobalVector(Gmtry%DataMngr,Cvt_zz_Global,ierr)

    CALL MPI_Comm_rank(MPI_COMM_WORLD,process,ierr)

    IF (process.EQ.0) THEN
        Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFileCvt))
        OPEN(u,FILE=TRIM(Route),STATUS='OLD',ACTION='READ')

        ! It gets the global size from the geometry data manager.
        CALL DMDAGetInfo(Gmtry%DataMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
        & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,ierr)

        DO i=1,widthG(1)*widthG(2)*widthG(3)
            READ(u,*)ValR(1),ValR(2),ValR(3)
            CALL VecSetValue(Cvt_xx_Global,i-1,ValR(1),INSERT_VALUES,ierr)
            CALL VecSetValue(Cvt_yy_Global,i-1,ValR(2),INSERT_VALUES,ierr)
            CALL VecSetValue(Cvt_zz_Global,i-1,ValR(3),INSERT_VALUES,ierr)
        END DO
        CLOSE(u)
    END IF

    CALL VecAssemblyBegin(Cvt_xx_Global,ierr)
    CALL VecAssemblyEnd(Cvt_xx_Global,ierr)
    CALL VecAssemblyBegin(Cvt_yy_Global,ierr)
    CALL VecAssemblyEnd(Cvt_yy_Global,ierr)
    CALL VecAssemblyBegin(Cvt_zz_Global,ierr)
    CALL VecAssemblyEnd(Cvt_zz_Global,ierr)

    EventName="GetConductivity"
    ViewName="Cvt_xx"
    CALL ViewConductivity(Cvt%xx,ViewName,EventName,ierr)
    CALL DMGlobalToLocalBegin(Gmtry%DataMngr,Cvt_xx_Global,INSERT_VALUES,Cvt%xx,ierr)
    CALL DMGlobalToLocalEnd(Gmtry%DataMngr,Cvt_xx_Global,INSERT_VALUES,Cvt%xx,ierr)
    CALL VecDestroy(Cvt_xx_Global,ierr)

    ViewName="Cvt_yy"
    CALL ViewConductivity(Cvt%yy,ViewName,EventName,ierr)
    CALL DMGlobalToLocalBegin(Gmtry%DataMngr,Cvt_yy_Global,INSERT_VALUES,Cvt%yy,ierr)
    CALL DMGlobalToLocalEnd(Gmtry%DataMngr,Cvt_yy_Global,INSERT_VALUES,Cvt%yy,ierr)
    CALL VecDestroy(Cvt_yy_Global,ierr)

    ViewName="Cvt_zz"
    CALL ViewConductivity(Cvt%zz,ViewName,EventName,ierr)
    CALL DMGlobalToLocalBegin(Gmtry%DataMngr,Cvt_zz_Global,INSERT_VALUES,Cvt%zz,ierr)
    CALL DMGlobalToLocalEnd(Gmtry%DataMngr,Cvt_zz_Global,INSERT_VALUES,Cvt%zz,ierr)
    CALL VecDestroy(Cvt_zz_Global,ierr)

    ! Saving a vector with zeros (If we put a conditionals on the stencil it can be avoided, but would be unclear)
    CALL DMCreateLocalVector(Gmtry%DataMngr,Cvt%xy,ierr)
    CALL VecZeroEntries(Cvt%xy,ierr)
    CALL VecAssemblyBegin(Cvt%xy,ierr)
    CALL VecAssemblyEnd(Cvt%xy,ierr)

    Cvt_xz_pointer => TargPetscVec(Cvt%xy)
    Cvt%xz = Cvt_xz_pointer
    Cvt_yz_pointer => TargPetscVec(Cvt%xy)
    Cvt%yz = Cvt_yz_pointer

!     IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetProrperties Stage] Conductivity Field was satisfactorily created\n",ierr)

END SUBROUTINE GetConductivity_3

SUBROUTINE GetConductivity_4(Gmtry,PptFld,Cvt,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,PropertiesField,ConductivityField
    USE ANISOFLOW_Interface, ONLY : GetVerbose,GetInputDir,GetInputFileCvt
    USE ANISOFLOW_View, ONLY : ViewConductivity
    
    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscdmda.h90>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    TYPE(Geometry),INTENT(IN)           :: Gmtry
    TYPE(PropertiesField),INTENT(IN)    :: PptFld
    TYPE(ConductivityField),INTENT(OUT) :: Cvt

    PetscInt                            :: u,i,widthG(3)
    PetscReal                           :: ValR(6)
    PetscMPIInt                         :: process
    CHARACTER(LEN=200)                  :: InputDir,InputFileCvt
    CHARACTER(LEN=200)                  :: Route,ViewName,EventName
    PetscBool                           :: Verbose
    Vec                                 :: Cvt_xx_Global,Cvt_yy_Global,Cvt_zz_Global
    Vec                                 :: Cvt_xy_Global,Cvt_xz_Global,Cvt_yz_Global

    PARAMETER(u=01)

    CALL GetVerbose(Verbose,ierr)

    CALL GetInputDir(InputDir,ierr)
    CALL GetInputFileCvt(InputFileCvt,ierr)

    Cvt%DefinedBy=3
!     IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetProrperties Event] Conductivity Field is stored by Block\n",ierr)

    CALL DMCreateLocalVector(Gmtry%DataMngr,Cvt%xx,ierr)
    CALL DMCreateLocalVector(Gmtry%DataMngr,Cvt%yy,ierr)
    CALL DMCreateLocalVector(Gmtry%DataMngr,Cvt%zz,ierr)
    CALL DMCreateLocalVector(Gmtry%DataMngr,Cvt%xy,ierr)
    CALL DMCreateLocalVector(Gmtry%DataMngr,Cvt%xz,ierr)
    CALL DMCreateLocalVector(Gmtry%DataMngr,Cvt%yz,ierr)

    CALL DMCreateGlobalVector(Gmtry%DataMngr,Cvt_xx_Global,ierr)
    CALL DMCreateGlobalVector(Gmtry%DataMngr,Cvt_yy_Global,ierr)
    CALL DMCreateGlobalVector(Gmtry%DataMngr,Cvt_zz_Global,ierr)
    CALL DMCreateGlobalVector(Gmtry%DataMngr,Cvt_xy_Global,ierr)
    CALL DMCreateGlobalVector(Gmtry%DataMngr,Cvt_xz_Global,ierr)
    CALL DMCreateGlobalVector(Gmtry%DataMngr,Cvt_yz_Global,ierr)

    CALL MPI_Comm_rank(MPI_COMM_WORLD,process,ierr)

    IF (process.EQ.0) THEN
        Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFileCvt))
        OPEN(u,FILE=TRIM(Route),STATUS='OLD',ACTION='READ')

        ! It gets the global size from the geometry data manager.
        CALL DMDAGetInfo(Gmtry%DataMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
        & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,ierr)

        DO i=1,widthG(1)*widthG(2)*widthG(3)
            READ(u,*)ValR(1),ValR(2),ValR(3),ValR(4),ValR(5),ValR(6)
            CALL VecSetValue(Cvt_xx_Global,i-1,ValR(1),INSERT_VALUES,ierr)
            CALL VecSetValue(Cvt_yy_Global,i-1,ValR(2),INSERT_VALUES,ierr)
            CALL VecSetValue(Cvt_zz_Global,i-1,ValR(3),INSERT_VALUES,ierr)
            CALL VecSetValue(Cvt_xy_Global,i-1,ValR(4),INSERT_VALUES,ierr)
            CALL VecSetValue(Cvt_xz_Global,i-1,ValR(5),INSERT_VALUES,ierr)
            CALL VecSetValue(Cvt_yz_Global,i-1,ValR(6),INSERT_VALUES,ierr)
        END DO
        CLOSE(u)
    END IF

    CALL VecAssemblyBegin(Cvt_xx_Global,ierr)
    CALL VecAssemblyEnd(Cvt_xx_Global,ierr)
    CALL VecAssemblyBegin(Cvt_yy_Global,ierr)
    CALL VecAssemblyEnd(Cvt_yy_Global,ierr)
    CALL VecAssemblyBegin(Cvt_zz_Global,ierr)
    CALL VecAssemblyEnd(Cvt_zz_Global,ierr)
    CALL VecAssemblyBegin(Cvt_xy_Global,ierr)
    CALL VecAssemblyEnd(Cvt_xy_Global,ierr)
    CALL VecAssemblyBegin(Cvt_xz_Global,ierr)
    CALL VecAssemblyEnd(Cvt_xz_Global,ierr)
    CALL VecAssemblyBegin(Cvt_yz_Global,ierr)
    CALL VecAssemblyEnd(Cvt_yz_Global,ierr)

    EventName="GetConductivity"
    ViewName="Cvt_xx"
    CALL ViewConductivity(Cvt%xx,ViewName,EventName,ierr)
    CALL DMGlobalToLocalBegin(Gmtry%DataMngr,Cvt_xx_Global,INSERT_VALUES,Cvt%xx,ierr)
    CALL DMGlobalToLocalEnd(Gmtry%DataMngr,Cvt_xx_Global,INSERT_VALUES,Cvt%xx,ierr)
    CALL VecDestroy(Cvt_xx_Global,ierr)

    ViewName="Cvt_yy"
    CALL ViewConductivity(Cvt%yy,ViewName,EventName,ierr)
    CALL DMGlobalToLocalBegin(Gmtry%DataMngr,Cvt_yy_Global,INSERT_VALUES,Cvt%yy,ierr)
    CALL DMGlobalToLocalEnd(Gmtry%DataMngr,Cvt_yy_Global,INSERT_VALUES,Cvt%yy,ierr)
    CALL VecDestroy(Cvt_yy_Global,ierr)

    ViewName="Cvt_zz"
    CALL ViewConductivity(Cvt%zz,ViewName,EventName,ierr)
    CALL DMGlobalToLocalBegin(Gmtry%DataMngr,Cvt_zz_Global,INSERT_VALUES,Cvt%zz,ierr)
    CALL DMGlobalToLocalEnd(Gmtry%DataMngr,Cvt_zz_Global,INSERT_VALUES,Cvt%zz,ierr)
    CALL VecDestroy(Cvt_zz_Global,ierr)

    ViewName="Cvt_xy"
    CALL ViewConductivity(Cvt%xy,ViewName,EventName,ierr)
    CALL DMGlobalToLocalBegin(Gmtry%DataMngr,Cvt_xy_Global,INSERT_VALUES,Cvt%xy,ierr)
    CALL DMGlobalToLocalEnd(Gmtry%DataMngr,Cvt_xy_Global,INSERT_VALUES,Cvt%xy,ierr)
    CALL VecDestroy(Cvt_xy_Global,ierr)

    ViewName="Cvt_xz"
    CALL ViewConductivity(Cvt%xz,ViewName,EventName,ierr)
    CALL DMGlobalToLocalBegin(Gmtry%DataMngr,Cvt_xz_Global,INSERT_VALUES,Cvt%xz,ierr)
    CALL DMGlobalToLocalEnd(Gmtry%DataMngr,Cvt_xz_Global,INSERT_VALUES,Cvt%xz,ierr)
    CALL VecDestroy(Cvt_xz_Global,ierr)

    ViewName="Cvt_yz"
    CALL ViewConductivity(Cvt%yz,ViewName,EventName,ierr)
    CALL DMGlobalToLocalBegin(Gmtry%DataMngr,Cvt_yz_Global,INSERT_VALUES,Cvt%yz,ierr)
    CALL DMGlobalToLocalEnd(Gmtry%DataMngr,Cvt_yz_Global,INSERT_VALUES,Cvt%yz,ierr)
    CALL VecDestroy(Cvt_yz_Global,ierr)

!     IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetProrperties Stage] Conductivity Field was satisfactorily created\n",ierr)

END SUBROUTINE GetConductivity_4

SUBROUTINE GetStorage(Gmtry,PptFld,Sto,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,PropertiesField,SpecificStorageField,InputTypeVar
    USE ANISOFLOW_Interface, ONLY : GetVerbose,GetInputType

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Geometry),INTENT(IN)               :: Gmtry
    TYPE(PropertiesField),INTENT(IN)        :: PptFld
    TYPE(SpecificStorageField),INTENT(OUT)  :: Sto

    TYPE(InputTypeVar)                  :: InputType
    PetscBool                           :: Verbose
    CHARACTER(LEN=200)                  :: EventName,ClassName
    PetscLogEvent                       :: Event
    PetscClassId                        :: ClassID
    PetscLogDouble                      :: EventFlops=0.d0

    ClassName="Property"
    CALL PetscClassIdRegister(ClassName,ClassID,ierr)
    EventName="GetSpecificStorage"
    CALL PetscLogEventRegister(EventName,ClassID,Event,ierr)
    CALL PetscLogEventBegin(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Inizialited\n",ierr)
    
    CALL GetInputType(InputType,ierr)

    IF (InputType%Sto.EQ.1) THEN
        CALL GetStorage_1(Gmtry,PptFld,Sto,ierr) ! By zone
        CALL StorageZoneToCell(Gmtry,Sto,ierr)
    ELSE IF (InputType%Sto.EQ.2) THEN
        CALL GetStorage_2(Gmtry,PptFld,Sto,ierr) ! By cell
    ELSE
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "[ERROR] Specific Storage InputType wrong\n",ierr)
        STOP
    END IF

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Finalized\n",ierr)
    
    CALL PetscLogFlops(EventFlops,ierr)
    CALL PetscLogEventEnd(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

END SUBROUTINE GetStorage

SUBROUTINE GetStoZoneID(Gmtry,PptFld,StoZoneID_Local,DefinedBy,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,PropertiesField,TargPetscVec
    USE ANISOFLOW_Interface, ONLY : GetVerbose,GetInputDir,GetInputFilePptByZone,GetInputFileStoByZone,GetHomogeneusPptFlg,GetHomogeneusCvtFlg
    USE ANISOFLOW_View, ONLY : ViewStorage
    USE ANISOFLOW_Geometry, ONLY : VecApplicationToPetsc

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscdmda.h90>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    TYPE(Geometry),INTENT(IN)           :: Gmtry
    TYPE(PropertiesField),INTENT(IN)    :: PptFld
    Vec,INTENT(OUT)                     :: StoZoneID_Local
    PetscInt,INTENT(OUT)                :: DefinedBy

    PetscInt                            :: widthG(3),u,ValI,i
    PetscMPIInt                         :: process
    PetscReal                           :: ValR,one=1.D0
    CHARACTER(LEN=200)                  :: InputFilePptByZone,InputFileStoByZone,Route,ViewName,InputDir,EventName
    PetscBool                           :: InputFilePptByZoneFlg,InputFileStoByZoneFlg,HomogeneusPptFlg,HomogeneusStoFlg,Verbose
    Vec,POINTER                         :: ZoneID_tmp
    Vec                                 :: StoZoneID_Global

    PARAMETER(u=01)

    CALL GetVerbose(Verbose,ierr)

    CALL GetInputDir(InputDir,ierr)
    CALL GetInputFilePptByZone(InputFilePptByZone,InputFilePptByZoneFlg,ierr)
    CALL GetInputFileStoByZone(InputFileStoByZone,InputFileStoByZoneFlg,ierr)

    CALL GetHomogeneusCvtFlg(HomogeneusStoFlg,ierr)
    CALL GetHomogeneusPptFlg(HomogeneusPptFlg,ierr)

    IF (InputFileStoByZoneFlg.AND.(.NOT.HomogeneusStoFlg)) THEN

        CALL DMCreateGlobalVector(Gmtry%DataMngr,StoZoneID_Global,ierr)
        CALL DMCreateLocalVector(Gmtry%DataMngr,StoZoneID_Local,ierr)

        CALL DMDAGetInfo(Gmtry%DataMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
            & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
            & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
            & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
            & PETSC_NULL_INTEGER,ierr)

        CALL MPI_Comm_rank(MPI_COMM_WORLD,process,ierr)
            IF (process.EQ.0) THEN
                Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFileStoByZone))
                OPEN(u,FILE=TRIM(Route),STATUS='OLD',ACTION='READ')
                DO i=1,widthG(1)*widthG(2)*widthG(3)
                    READ(u,*)ValI
                    IF (ValI.LE.0) THEN
                        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,             &
                            & "[ERROR] Input_file_sto_by_zones entry only can contain"  &
                            & // " natural numbers\n",ierr)
                        STOP
                    END IF
                    ValR=REAL(ValI)
                    CALL VecSetValue(StoZoneID_Local,i-1,ValR,INSERT_VALUES,ierr)
                END DO
                CLOSE(u)
            END IF

        CALL VecAssemblyBegin(StoZoneID_Local,ierr)
        CALL VecAssemblyEnd(StoZoneID_Local,ierr)

        CALL VecApplicationToPetsc(Gmtry%DataMngr,StoZoneID_Local,ierr)

        ViewName="StoZoneID"
        EventName="GetStorage"
        CALL ViewStorage(StoZoneID_Global,ViewName,EventName,ierr)

        CALL DMGlobalToLocalBegin(Gmtry%DataMngr,StoZoneID_Global,INSERT_VALUES,StoZoneID_Local,ierr)
        CALL DMGlobalToLocalEnd(Gmtry%DataMngr,StoZoneID_Global,INSERT_VALUES,StoZoneID_Local,ierr)
        CALL VecDestroy(StoZoneID_Global,ierr)
        DefinedBy=1
    ELSEIF (HomogeneusStoFlg) THEN

        CALL DMCreateGlobalVector(Gmtry%DataMngr,StoZoneID_Global,ierr)
        CALL DMCreateLocalVector(Gmtry%DataMngr,StoZoneID_Local,ierr)

        CALL DMDAGetInfo(Gmtry%DataMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
            & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
            & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
            & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
            & PETSC_NULL_INTEGER,ierr)

        CALL MPI_Comm_rank(MPI_COMM_WORLD,process,ierr)

        CALL VecSet(StoZoneID_Global,one,ierr)

        CALL VecAssemblyBegin(StoZoneID_Global,ierr)
        CALL VecAssemblyEnd(StoZoneID_Global,ierr)

        CALL VecApplicationToPetsc(Gmtry%DataMngr,StoZoneID_Global,ierr)

        ViewName="StoZoneID"
        EventName="GetStoZoneID"
        CALL ViewStorage(StoZoneID_Global,ViewName,EventName,ierr)

        CALL DMGlobalToLocalBegin(Gmtry%DataMngr,StoZoneID_Global,INSERT_VALUES,StoZoneID_Local,ierr)
        CALL DMGlobalToLocalEnd(Gmtry%DataMngr,StoZoneID_Global,INSERT_VALUES,StoZoneID_Local,ierr)
        CALL VecDestroy(StoZoneID_Global,ierr)

        DefinedBy=1
        IF (InputFileStoByZoneFlg) THEN 
            IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetProrperties Event] WARNING: -Input_file_sto_by_zones and -Homogeneuos_sto used at the same time; specific storages was set homogeneuos in the whole domain.\n",ierr)
        END IF
    ELSEIF (InputFilePptByZoneFlg.OR.HomogeneusPptFlg) THEN
        ZoneID_tmp => TargPetscVec(PptFld%ZoneID)
        StoZoneID_Local = ZoneID_tmp
        DefinedBy=2
    ELSE
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,             &
            & "[ERROR] Input_type_sto must use Input_file_sto_by_zones or Input_file_ppt_by_zones"  &
            & // " containing the zonification of specific storage.\n",ierr)
        STOP
    END IF

END SUBROUTINE GetStoZoneID

SUBROUTINE GetStorage_1(Gmtry,PptFld,Sto,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,SpecificStorageField,PropertiesField,TargetFullTensor
    USE ANISOFLOW_Interface, ONLY : GetVerbose,GetInputDir,GetInputFileSto
    USE ANISOFLOW_Geometry, ONLY : VecApplicationToPetsc
    
    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscdmda.h90>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Geometry),INTENT(IN)               :: Gmtry
    TYPE(PropertiesField),INTENT(IN)        :: PptFld
    TYPE(SpecificStorageField),INTENT(OUT)  :: Sto

    PetscInt                            :: u,i,j,StoLen
    PetscReal                           :: ValR
    PetscMPIInt                         :: process
    CHARACTER(LEN=200)                  :: InputDir,InputFileSto!,ViewName!,EventName
    CHARACTER(LEN=200)                  :: Route
    CHARACTER(LEN=16)                   :: StoKind
    PetscBool                           :: Verbose
    PetscReal,ALLOCATABLE               :: Zone_tmp(:)

    PARAMETER(u=01)

    CALL GetVerbose(Verbose,ierr)

    CALL GetStoZoneID(Gmtry,PptFld,Sto%ZoneID,Sto%DefinedBy,ierr)
    CALL GetInputDir(InputDir,ierr)
    CALL GetInputFileSto(InputFileSto,ierr)

    CALL MPI_Comm_rank(MPI_COMM_WORLD,process,ierr)

!     IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetStorage Event] Storage Field is stored as Zones by Block\n",ierr)

    CALL VecMax(Sto%ZoneID,PETSC_NULL_INTEGER,ValR,ierr)
    StoLen=INT(ValR)

    ALLOCATE(Zone_tmp(StoLen))

    IF (process.EQ.0) THEN
        Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFileSto))
        OPEN(u,FILE=TRIM(Route),STATUS='OLD',ACTION='READ')

        DO i=1,StoLen
            DO j=1,6
                READ(u,*)
            END DO
            READ(u,'(A16)')StoKind
            IF (StoKind.EQ."specific storage") THEN 
                READ(u,*)Zone_tmp(i)
            ELSE
                CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,             &
                    & "[ERROR] File of Specific Storage properties is invalid\n"&
                    & ,ierr)
                STOP            
            END IF
            DO j=1,21
                READ(u,*)
            END DO
        END DO
        CLOSE(u)
    END IF


    CALL MPI_Bcast(Zone_tmp(:),StoLen,MPI_DOUBLE,0,PETSC_COMM_WORLD,ierr)

    CALL VecCreateSeq(PETSC_COMM_SELF,StoLen,Sto%Zone,ierr)

    CALL VecSetValues(Sto%Zone,StoLen,(/(i,i=0,StoLen-1)/),Zone_tmp,INSERT_VALUES,ierr)

    CALL VecAssemblyBegin(Sto%Zone,ierr)
    CALL VecAssemblyEnd(Sto%Zone,ierr)

!     IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetStorage Event] Storage Field was satisfactorily created\n",ierr)

END SUBROUTINE GetStorage_1

SUBROUTINE StorageZoneToCell(Gmtry,Sto,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,SpecificStorageField
    USE ANISOFLOW_View, ONLY : ViewStorage

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscdmda.h90>

    PetscErrorCode,INTENT(INOUT)                :: ierr
    TYPE(Geometry),INTENT(IN)                   :: Gmtry
    TYPE(SpecificStorageField),INTENT(INOUT)    :: Sto

    CHARACTER(LEN=200)      :: ViewName,EventName
    PetscInt                :: widthL(3),widthG(3),corn(3),ValI,i,j,k
    PetscReal,POINTER       :: TmpStoZone(:),TmpStoZoneID3D(:,:,:),TmpStoZoneID2D(:,:),TmpStoCell3D(:,:,:),TmpStoCell2D(:,:)

    CALL DMDAGetInfo(Gmtry%DataMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
        & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,ierr)

    IF (.NOT.((Sto%DefinedBy.EQ.1).OR.(Sto%DefinedBy.EQ.2))) THEN
        PRINT*,"ERROR: To use StorageZoneToCell subroutine the specific storage has to be definde by zone."
    END IF

    CALL DMDAGetCorners(Gmtry%DataMngr,corn(1),corn(2),corn(3),widthL(1),widthL(2),widthL(3),ierr)
    CALL DMCreateGlobalVector(Gmtry%DataMngr,Sto%Cell,ierr)

    IF (widthG(3).NE.1) THEN
        CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,Sto%ZoneID,TmpStoZoneID3D,ierr)
        CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,Sto%Cell,TmpStoCell3D,ierr)
    ELSE
        CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,Sto%ZoneID,TmpStoZoneID2D,ierr)
        CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,Sto%Cell,TmpStoCell2D,ierr)
    END IF       

    CALL VecGetArrayReadF90(Sto%Zone,TmpStoZone,ierr)

    DO k=corn(3),corn(3)+widthL(3)-1
        DO j=corn(2),corn(2)+widthL(2)-1
            DO i=corn(1),corn(1)+widthL(1)-1
                IF (widthG(3).NE.1) THEN
                    ValI=INT(TmpStoZoneID3D(i,j,k))
                    TmpStoCell3D(i,j,k)=TmpStoZone(ValI)
                ELSE
                    ValI=INT(TmpStoZoneID2D(i,j))
                    TmpStoCell2D(i,j)=TmpStoZone(ValI)
                END IF
            END DO
        END DO
    END DO

    CALL VecRestoreArrayReadF90(Sto%Zone,TmpStoZone,ierr)
    IF (widthG(3).NE.1) THEN
        CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,Sto%ZoneID,TmpStoZoneID3D,ierr)
        CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,Sto%Cell,TmpStoCell3D,ierr)
    ELSE
        CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,Sto%ZoneID,TmpStoZoneID2D,ierr)
        CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,Sto%Cell,TmpStoCell2D,ierr)
    END IF
    
    CALL VecDestroy(Sto%Zone,ierr)

    IF (Sto%DefinedBy.EQ.1) THEN
        CALL VecDestroy(Sto%ZoneID,ierr)
        Sto%DefinedBy=0
    ELSE IF (Sto%DefinedBy.EQ.2) THEN
        !NULLIFY(Sto%ZoneID) ! Because I used a very tricky way to use it as a pointer, I can't nullify it as a usual pointer but it's supposed it doesn't will get in troubles the code.
        Sto%DefinedBy=0
    END IF

    Sto%DefinedBy=3

    ViewName="Sto"
    EventName="StorageZoneToCell"
    CALL ViewStorage(Sto%Cell,ViewName,EventName,ierr)

END SUBROUTINE StorageZoneToCell

SUBROUTINE GetStorage_2(Gmtry,PptFld,Sto,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,PropertiesField,SpecificStorageField
    USE ANISOFLOW_Interface, ONLY : GetVerbose,GetInputDir,GetInputFileSto
    USE ANISOFLOW_View, ONLY : ViewStorage
    
    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscdmda.h90>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Geometry),INTENT(IN)               :: Gmtry
    TYPE(PropertiesField),INTENT(IN)        :: PptFld
    TYPE(SpecificStorageField),INTENT(OUT)  :: Sto

    PetscInt                            :: u,i,widthG(3)
    PetscReal                           :: ValR
    PetscMPIInt                         :: process
    CHARACTER(LEN=200)                  :: InputDir,InputFileSto
    CHARACTER(LEN=200)                  :: Route,ViewName,EventName
    PetscBool                           :: Verbose

    PARAMETER(u=01)

    CALL GetVerbose(Verbose,ierr)

    CALL GetInputDir(InputDir,ierr)
    CALL GetInputFileSto(InputFileSto,ierr)

    Sto%DefinedBy=3
!     IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetProrperties Event] Storage Field is stored by Block\n",ierr)

    CALL DMCreateGlobalVector(Gmtry%DataMngr,Sto%Cell,ierr)

    CALL MPI_Comm_rank(MPI_COMM_WORLD,process,ierr)

    IF (process.EQ.0) THEN
        Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFileSto))
        OPEN(u,FILE=TRIM(Route),STATUS='OLD',ACTION='READ')

        ! It gets the global size from the geometry data manager.
        CALL DMDAGetInfo(Gmtry%DataMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
        & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,ierr)

        DO i=1,widthG(1)*widthG(2)*widthG(3)
            READ(u,*)ValR
            CALL VecSetValue(Sto%Cell,i-1,ValR,INSERT_VALUES,ierr)
        END DO
        CLOSE(u)
    END IF

    CALL VecAssemblyBegin(Sto%Cell,ierr)
    CALL VecAssemblyEnd(Sto%Cell,ierr)

    ViewName="Sto"
    EventName="GetStorage"
    CALL ViewStorage(Sto%Cell,ViewName,EventName,ierr)

!     IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetProrperties Stage] Storage Field was satisfactorily created\n",ierr)

END SUBROUTINE GetStorage_2


SUBROUTINE GetLocalProperty(Gmtry,PptFld,Ppt,i,j,k,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,PropertiesField,Property
    USE ANISOFLOW_Geometry, ONLY : GetLocalTopology

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    TYPE(Geometry),INTENT(IN)           :: Gmtry
    TYPE(PropertiesField),INTENT(IN)    :: PptFld
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

    USE ANISOFLOW_Types, ONLY : Geometry,PropertiesField,Property,Tensor,TargetFullTensor
    USE ANISOFLOW_Operators

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscdmda.h90>

    TYPE(Geometry),INTENT(IN)           :: Gmtry
    TYPE(PropertiesField),INTENT(IN)    :: PptFld
    TYPE(Property),INTENT(INOUT)        :: Ppt
    PetscErrorCode,INTENT(INOUT)        :: ierr

    PetscReal,POINTER                   :: TmpCvt3D(:,:,:),TmpCvt2D(:,:)
    TYPE(Tensor)                        :: TensorZero
    PetscInt                            :: CenterID,widthG(3),i,j,k
    PetscInt,ALLOCATABLE                :: ValI(:)

    i=Ppt%Pstn%i
    j=Ppt%Pstn%j
    k=Ppt%Pstn%k

    ! It gets the global size from the geometry data manager.
    CALL DMDAGetInfo(Gmtry%DataMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
        & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,ierr)

    IF (Ppt%StnclType.EQ.1) THEN          ! Stencil type star
        ALLOCATE(Ppt%Cvt(7),ValI(7))
        CenterID=4
    ELSEIF (Ppt%StnclType.EQ.2) THEN      ! Stencil type box (19)
        ALLOCATE(Ppt%Cvt(19),ValI(19))
        CenterID=10
    ELSE 
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,              &
            & "[ERROR] Property topology wasn't well defined\n"      &
            & ,ierr)
        STOP   
    END IF
    
    TensorZero%xx=0.D0
    TensorZero%yy=0.D0
    TensorZero%zz=0.D0
    TensorZero%xy=0.D0
    TensorZero%xz=0.D0
    TensorZero%yz=0.D0
     CALL TargetFullTensor(TensorZero)

    Ppt%Cvt(:)=TensorZero

    IF (Ppt%StnclType.EQ.1) THEN          ! Stencil type star
        IF ((PptFld%Cvt%DefinedBy.EQ.1).OR.(PptFld%Cvt%DefinedBy.EQ.2)) THEN 
            ! Defined by zones, so have to get the information from Zone(ZoneID(i))
            IF ((Ppt%StnclTplgy(CenterID).EQ.1).OR.(Ppt%StnclTplgy(CenterID).EQ.3).OR.(Ppt%StnclTplgy(CenterID).EQ.4)) THEN 
                ! Only get properties for active blocks
                IF (widthG(3).NE.1) THEN ! 3D mesh
                    
                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%ZoneID,TmpCvt3D,ierr)
                    ValI(1)=INT(TmpCvt3D(i  ,j  ,k-1))
                    ValI(2)=INT(TmpCvt3D(i  ,j-1,k  ))
                    ValI(3)=INT(TmpCvt3D(i-1,j  ,k  ))
                    ValI(4)=INT(TmpCvt3D(i  ,j  ,k  ))
                    ValI(5)=INT(TmpCvt3D(i+1,j  ,k  ))
                    ValI(6)=INT(TmpCvt3D(i  ,j+1,k  ))
                    ValI(7)=INT(TmpCvt3D(i  ,j  ,k+1))
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%ZoneID,TmpCvt3D,ierr)
                    Ppt%Cvt(:)=PptFld%Cvt%Zone(ValI(:))
                ELSE ! 2D mesh
                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%ZoneID,TmpCvt2D,ierr)
                    ValI(2)=INT(TmpCvt2D(i  ,j-1))
                    ValI(3)=INT(TmpCvt2D(i-1,j  ))
                    ValI(4)=INT(TmpCvt2D(i  ,j  ))
                    ValI(5)=INT(TmpCvt2D(i+1,j  ))
                    ValI(6)=INT(TmpCvt2D(i  ,j+1))
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%ZoneID,TmpCvt2D,ierr)
                    Ppt%Cvt(2:6)=PptFld%Cvt%Zone(ValI(2:6))
                END IF                
            END IF
        ELSE IF (Pptfld%Cvt%DefinedBy.EQ.3) THEN
            IF ((Ppt%StnclTplgy(CenterID).EQ.1).OR.(Ppt%StnclTplgy(CenterID).EQ.3).OR.(Ppt%StnclTplgy(CenterID).EQ.4)) THEN

                IF (widthG(3).NE.1) THEN  ! 3D mesh
                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%xx,TmpCvt3D,ierr)
                    Ppt%Cvt(1)%xx=TmpCvt3D(i  ,j  ,k-1)
                    Ppt%Cvt(2)%xx=TmpCvt3D(i  ,j-1,k  )
                    Ppt%Cvt(3)%xx=TmpCvt3D(i-1,j  ,k  )
                    Ppt%Cvt(4)%xx=TmpCvt3D(i  ,j  ,k  )
                    Ppt%Cvt(5)%xx=TmpCvt3D(i+1,j  ,k  )
                    Ppt%Cvt(6)%xx=TmpCvt3D(i  ,j+1,k  )
                    Ppt%Cvt(7)%xx=TmpCvt3D(i  ,j  ,k+1)
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%xx,TmpCvt3D,ierr)

                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%yy,TmpCvt3D,ierr)
                    Ppt%Cvt(1)%yy=TmpCvt3D(i  ,j  ,k-1)
                    Ppt%Cvt(2)%yy=TmpCvt3D(i  ,j-1,k  )
                    Ppt%Cvt(3)%yy=TmpCvt3D(i-1,j  ,k  )
                    Ppt%Cvt(4)%yy=TmpCvt3D(i  ,j  ,k  )
                    Ppt%Cvt(5)%yy=TmpCvt3D(i+1,j  ,k  )
                    Ppt%Cvt(6)%yy=TmpCvt3D(i  ,j+1,k  )
                    Ppt%Cvt(7)%yy=TmpCvt3D(i  ,j  ,k+1)
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%yy,TmpCvt3D,ierr)

                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%zz,TmpCvt3D,ierr)
                    Ppt%Cvt(1)%zz=TmpCvt3D(i  ,j  ,k-1)
                    Ppt%Cvt(2)%zz=TmpCvt3D(i  ,j-1,k  )
                    Ppt%Cvt(3)%zz=TmpCvt3D(i-1,j  ,k  )
                    Ppt%Cvt(4)%zz=TmpCvt3D(i  ,j  ,k  )
                    Ppt%Cvt(5)%zz=TmpCvt3D(i+1,j  ,k  )
                    Ppt%Cvt(6)%zz=TmpCvt3D(i  ,j+1,k  )
                    Ppt%Cvt(7)%zz=TmpCvt3D(i  ,j  ,k+1)
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%zz,TmpCvt3D,ierr)

                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%xy,TmpCvt3D,ierr)
                    Ppt%Cvt(1)%xy=TmpCvt3D(i  ,j  ,k-1)
                    Ppt%Cvt(2)%xy=TmpCvt3D(i  ,j-1,k  )
                    Ppt%Cvt(3)%xy=TmpCvt3D(i-1,j  ,k  )
                    Ppt%Cvt(4)%xy=TmpCvt3D(i  ,j  ,k  )
                    Ppt%Cvt(5)%xy=TmpCvt3D(i+1,j  ,k  )
                    Ppt%Cvt(6)%xy=TmpCvt3D(i  ,j+1,k  )
                    Ppt%Cvt(7)%xy=TmpCvt3D(i  ,j  ,k+1)
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%xy,TmpCvt3D,ierr)

                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%xz,TmpCvt3D,ierr)
                    Ppt%Cvt(1)%xz=TmpCvt3D(i  ,j  ,k-1)
                    Ppt%Cvt(2)%xz=TmpCvt3D(i  ,j-1,k  )
                    Ppt%Cvt(3)%xz=TmpCvt3D(i-1,j  ,k  )
                    Ppt%Cvt(4)%xz=TmpCvt3D(i  ,j  ,k  )
                    Ppt%Cvt(5)%xz=TmpCvt3D(i+1,j  ,k  )
                    Ppt%Cvt(6)%xz=TmpCvt3D(i  ,j+1,k  )
                    Ppt%Cvt(7)%xz=TmpCvt3D(i  ,j  ,k+1)
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%xz,TmpCvt3D,ierr)

                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%yz,TmpCvt3D,ierr)
                    Ppt%Cvt(1)%yz=TmpCvt3D(i  ,j  ,k-1)
                    Ppt%Cvt(2)%yz=TmpCvt3D(i  ,j-1,k  )
                    Ppt%Cvt(3)%yz=TmpCvt3D(i-1,j  ,k  )
                    Ppt%Cvt(4)%yz=TmpCvt3D(i  ,j  ,k  )
                    Ppt%Cvt(5)%yz=TmpCvt3D(i+1,j  ,k  )
                    Ppt%Cvt(6)%yz=TmpCvt3D(i  ,j+1,k  )
                    Ppt%Cvt(7)%yz=TmpCvt3D(i  ,j  ,k+1)
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%yz,TmpCvt3D,ierr)

                    CALL TargetFullTensor(Ppt%Cvt(1))
                    CALL TargetFullTensor(Ppt%Cvt(2))
                    CALL TargetFullTensor(Ppt%Cvt(3))
                    CALL TargetFullTensor(Ppt%Cvt(4))
                    CALL TargetFullTensor(Ppt%Cvt(5))
                    CALL TargetFullTensor(Ppt%Cvt(6))
                    CALL TargetFullTensor(Ppt%Cvt(7))

                ELSE ! 2D mesh
                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%xx,TmpCvt2D,ierr)
                    Ppt%Cvt(2)%xx=TmpCvt2D(i  ,j-1)
                    Ppt%Cvt(3)%xx=TmpCvt2D(i-1,j  )
                    Ppt%Cvt(4)%xx=TmpCvt2D(i  ,j  )
                    Ppt%Cvt(5)%xx=TmpCvt2D(i+1,j  )
                    Ppt%Cvt(6)%xx=TmpCvt2D(i  ,j+1)
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%xx,TmpCvt2D,ierr)

                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%yy,TmpCvt2D,ierr)
                    Ppt%Cvt(2)%yy=TmpCvt2D(i  ,j-1)
                    Ppt%Cvt(3)%yy=TmpCvt2D(i-1,j  )
                    Ppt%Cvt(4)%yy=TmpCvt2D(i  ,j  )
                    Ppt%Cvt(5)%yy=TmpCvt2D(i+1,j  )
                    Ppt%Cvt(6)%yy=TmpCvt2D(i  ,j+1)
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%yy,TmpCvt2D,ierr)

                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%zz,TmpCvt2D,ierr)
                    Ppt%Cvt(2)%zz=TmpCvt2D(i  ,j-1)
                    Ppt%Cvt(3)%zz=TmpCvt2D(i-1,j  )
                    Ppt%Cvt(4)%zz=TmpCvt2D(i  ,j  )
                    Ppt%Cvt(5)%zz=TmpCvt2D(i+1,j  )
                    Ppt%Cvt(6)%zz=TmpCvt2D(i  ,j+1)
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%zz,TmpCvt2D,ierr)

                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%xy,TmpCvt2D,ierr)
                    Ppt%Cvt(2)%xy=TmpCvt2D(i  ,j-1)
                    Ppt%Cvt(3)%xy=TmpCvt2D(i-1,j  )
                    Ppt%Cvt(4)%xy=TmpCvt2D(i  ,j  )
                    Ppt%Cvt(5)%xy=TmpCvt2D(i+1,j  )
                    Ppt%Cvt(6)%xy=TmpCvt2D(i  ,j+1)
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%xy,TmpCvt2D,ierr)

                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%xz,TmpCvt2D,ierr)
                    Ppt%Cvt(2)%xz=TmpCvt2D(i  ,j-1)
                    Ppt%Cvt(3)%xz=TmpCvt2D(i-1,j  )
                    Ppt%Cvt(4)%xz=TmpCvt2D(i  ,j  )
                    Ppt%Cvt(5)%xz=TmpCvt2D(i+1,j  )
                    Ppt%Cvt(6)%xz=TmpCvt2D(i  ,j+1)
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%xz,TmpCvt2D,ierr)

                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%yz,TmpCvt2D,ierr)
                    Ppt%Cvt(2)%yz=TmpCvt2D(i  ,j-1)
                    Ppt%Cvt(3)%yz=TmpCvt2D(i-1,j  )
                    Ppt%Cvt(4)%yz=TmpCvt2D(i  ,j  )
                    Ppt%Cvt(5)%yz=TmpCvt2D(i+1,j  )
                    Ppt%Cvt(6)%yz=TmpCvt2D(i  ,j+1)
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%yz,TmpCvt2D,ierr)

                    CALL TargetFullTensor(Ppt%Cvt(2))
                    CALL TargetFullTensor(Ppt%Cvt(3))
                    CALL TargetFullTensor(Ppt%Cvt(4))
                    CALL TargetFullTensor(Ppt%Cvt(5))
                    CALL TargetFullTensor(Ppt%Cvt(6))

                END IF      
            ELSE
                PRINT*,"ERROR: The properties has to be defined by cell or by zone"
            END IF  
        END IF
        IF (k.EQ.0) Ppt%Cvt(1)=TensorZero
        IF (j.EQ.0) Ppt%Cvt(2)=TensorZero
        IF (i.EQ.0) Ppt%Cvt(3)=TensorZero
        IF (i.EQ.(widthG(1)-1)) Ppt%Cvt(5)=TensorZero
        IF (j.EQ.(widthG(2)-1)) Ppt%Cvt(6)=TensorZero
        IF (k.EQ.(widthG(3)-1)) Ppt%Cvt(7)=TensorZero

    ELSEIF (Ppt%StnclType.EQ.2) THEN          ! Stencil type box (19)
        IF ((PptFld%Cvt%DefinedBy.EQ.1).OR.(PptFld%Cvt%DefinedBy.EQ.2)) THEN 
            ! Defined by zones, so have to get the information from Zone(ZoneID(i))
            IF ((Ppt%StnclTplgy(CenterID).EQ.1).OR.(Ppt%StnclTplgy(CenterID).EQ.3).OR.(Ppt%StnclTplgy(CenterID).EQ.4)) THEN 
                ! Only get properties for active blocks
                IF (widthG(3).NE.1) THEN ! 3D mesh
                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%ZoneID,TmpCvt3D,ierr)
                    ValI(1)= INT(TmpCvt3D(i  ,j-1,k-1))
                    ValI(2)= INT(TmpCvt3D(i-1,j  ,k-1))
                    ValI(3)= INT(TmpCvt3D(i  ,j  ,k-1))
                    ValI(4)= INT(TmpCvt3D(i+1,j  ,k-1))
                    ValI(5)= INT(TmpCvt3D(i  ,j+1,k-1))
                    ValI(6)= INT(TmpCvt3D(i-1,j-1,k  ))
                    ValI(7)= INT(TmpCvt3D(i  ,j-1,k  ))
                    ValI(8)= INT(TmpCvt3D(i+1,j-1,k  ))
                    ValI(9)= INT(TmpCvt3D(i-1,j  ,k  ))
                    ValI(10)=INT(TmpCvt3D(i  ,j  ,k  ))
                    ValI(11)=INT(TmpCvt3D(i+1,j  ,k  ))
                    ValI(12)=INT(TmpCvt3D(i-1,j+1,k  ))
                    ValI(13)=INT(TmpCvt3D(i  ,j+1,k  ))
                    ValI(14)=INT(TmpCvt3D(i+1,j+1,k  ))
                    ValI(15)=INT(TmpCvt3D(i  ,j-1,k+1))
                    ValI(16)=INT(TmpCvt3D(i-1,j  ,k+1))
                    ValI(17)=INT(TmpCvt3D(i  ,j  ,k+1))
                    ValI(18)=INT(TmpCvt3D(i+1,j  ,k+1))
                    ValI(19)=INT(TmpCvt3D(i  ,j+1,k+1))
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%ZoneID,TmpCvt3D,ierr)
                    Ppt%Cvt(:)=PptFld%Cvt%Zone(ValI(:))

                ELSE ! 2D mesh
                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%ZoneID,TmpCvt2D,ierr)
                    ValI(6)= INT(TmpCvt2D(i-1,j-1))
                    ValI(7)= INT(TmpCvt2D(i  ,j-1))
                    ValI(8)= INT(TmpCvt2D(i+1,j-1))
                    ValI(9)= INT(TmpCvt2D(i-1,j  ))
                    ValI(10)=INT(TmpCvt2D(i  ,j  ))
                    ValI(11)=INT(TmpCvt2D(i+1,j  ))
                    ValI(12)=INT(TmpCvt2D(i-1,j+1))
                    ValI(13)=INT(TmpCvt2D(i  ,j+1))
                    ValI(14)=INT(TmpCvt2D(i+1,j+1))
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%ZoneID,TmpCvt2D,ierr)
                    Ppt%Cvt(6:14)=PptFld%Cvt%Zone(ValI(6:14))
                END IF                
            END IF
        ELSE IF (Pptfld%Cvt%DefinedBy.EQ.3) THEN
            IF ((Ppt%StnclTplgy(CenterID).EQ.1).OR.(Ppt%StnclTplgy(CenterID).EQ.3).OR.(Ppt%StnclTplgy(CenterID).EQ.4)) THEN

                IF (widthG(3).NE.1) THEN  ! 3D mesh
                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%xx,TmpCvt3D,ierr)
                    Ppt%Cvt(1)%xx= TmpCvt3D(i  ,j-1,k-1)
                    Ppt%Cvt(2)%xx= TmpCvt3D(i-1,j  ,k-1)
                    Ppt%Cvt(3)%xx= TmpCvt3D(i  ,j  ,k-1)
                    Ppt%Cvt(4)%xx= TmpCvt3D(i+1,j  ,k-1)
                    Ppt%Cvt(5)%xx= TmpCvt3D(i  ,j+1,k-1)
                    Ppt%Cvt(6)%xx= TmpCvt3D(i-1,j-1,k  )
                    Ppt%Cvt(7)%xx= TmpCvt3D(i  ,j-1,k  )
                    Ppt%Cvt(8)%xx= TmpCvt3D(i+1,j-1,k  )
                    Ppt%Cvt(9)%xx= TmpCvt3D(i-1,j  ,k  )
                    Ppt%Cvt(10)%xx=TmpCvt3D(i  ,j  ,k  )
                    Ppt%Cvt(11)%xx=TmpCvt3D(i+1,j  ,k  )
                    Ppt%Cvt(12)%xx=TmpCvt3D(i-1,j+1,k  )
                    Ppt%Cvt(13)%xx=TmpCvt3D(i  ,j+1,k  )
                    Ppt%Cvt(14)%xx=TmpCvt3D(i+1,j+1,k  )
                    Ppt%Cvt(15)%xx=TmpCvt3D(i  ,j-1,k+1)
                    Ppt%Cvt(16)%xx=TmpCvt3D(i-1,j  ,k+1)
                    Ppt%Cvt(17)%xx=TmpCvt3D(i  ,j  ,k+1)
                    Ppt%Cvt(18)%xx=TmpCvt3D(i+1,j  ,k+1)
                    Ppt%Cvt(19)%xx=TmpCvt3D(i  ,j+1,k+1)
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%xx,TmpCvt3D,ierr)

                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%yy,TmpCvt3D,ierr)
                    Ppt%Cvt(1)%yy= TmpCvt3D(i  ,j-1,k-1)
                    Ppt%Cvt(2)%yy= TmpCvt3D(i-1,j  ,k-1)
                    Ppt%Cvt(3)%yy= TmpCvt3D(i  ,j  ,k-1)
                    Ppt%Cvt(4)%yy= TmpCvt3D(i+1,j  ,k-1)
                    Ppt%Cvt(5)%yy= TmpCvt3D(i  ,j+1,k-1)
                    Ppt%Cvt(6)%yy= TmpCvt3D(i-1,j-1,k  )
                    Ppt%Cvt(7)%yy= TmpCvt3D(i  ,j-1,k  )
                    Ppt%Cvt(8)%yy= TmpCvt3D(i+1,j-1,k  )
                    Ppt%Cvt(9)%yy= TmpCvt3D(i-1,j  ,k  )
                    Ppt%Cvt(10)%yy=TmpCvt3D(i  ,j  ,k  )
                    Ppt%Cvt(11)%yy=TmpCvt3D(i+1,j  ,k  )
                    Ppt%Cvt(12)%yy=TmpCvt3D(i-1,j+1,k  )
                    Ppt%Cvt(13)%yy=TmpCvt3D(i  ,j+1,k  )
                    Ppt%Cvt(14)%yy=TmpCvt3D(i+1,j+1,k  )
                    Ppt%Cvt(15)%yy=TmpCvt3D(i  ,j-1,k+1)
                    Ppt%Cvt(16)%yy=TmpCvt3D(i-1,j  ,k+1)
                    Ppt%Cvt(17)%yy=TmpCvt3D(i  ,j  ,k+1)
                    Ppt%Cvt(18)%yy=TmpCvt3D(i+1,j  ,k+1)
                    Ppt%Cvt(19)%yy=TmpCvt3D(i  ,j+1,k+1)
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%yy,TmpCvt3D,ierr)

                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%zz,TmpCvt3D,ierr)
                    Ppt%Cvt(1)%zz= TmpCvt3D(i  ,j-1,k-1)
                    Ppt%Cvt(2)%zz= TmpCvt3D(i-1,j  ,k-1)
                    Ppt%Cvt(3)%zz= TmpCvt3D(i  ,j  ,k-1)
                    Ppt%Cvt(4)%zz= TmpCvt3D(i+1,j  ,k-1)
                    Ppt%Cvt(5)%zz= TmpCvt3D(i  ,j+1,k-1)
                    Ppt%Cvt(6)%zz= TmpCvt3D(i-1,j-1,k  )
                    Ppt%Cvt(7)%zz= TmpCvt3D(i  ,j-1,k  )
                    Ppt%Cvt(8)%zz= TmpCvt3D(i+1,j-1,k  )
                    Ppt%Cvt(9)%zz= TmpCvt3D(i-1,j  ,k  )
                    Ppt%Cvt(10)%zz=TmpCvt3D(i  ,j  ,k  )
                    Ppt%Cvt(11)%zz=TmpCvt3D(i+1,j  ,k  )
                    Ppt%Cvt(12)%zz=TmpCvt3D(i-1,j+1,k  )
                    Ppt%Cvt(13)%zz=TmpCvt3D(i  ,j+1,k  )
                    Ppt%Cvt(14)%zz=TmpCvt3D(i+1,j+1,k  )
                    Ppt%Cvt(15)%zz=TmpCvt3D(i  ,j-1,k+1)
                    Ppt%Cvt(16)%zz=TmpCvt3D(i-1,j  ,k+1)
                    Ppt%Cvt(17)%zz=TmpCvt3D(i  ,j  ,k+1)
                    Ppt%Cvt(18)%zz=TmpCvt3D(i+1,j  ,k+1)
                    Ppt%Cvt(19)%zz=TmpCvt3D(i  ,j+1,k+1)
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%zz,TmpCvt3D,ierr)

                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%xy,TmpCvt3D,ierr)
                    Ppt%Cvt(1)%xy= TmpCvt3D(i  ,j-1,k-1)
                    Ppt%Cvt(2)%xy= TmpCvt3D(i-1,j  ,k-1)
                    Ppt%Cvt(3)%xy= TmpCvt3D(i  ,j  ,k-1)
                    Ppt%Cvt(4)%xy= TmpCvt3D(i+1,j  ,k-1)
                    Ppt%Cvt(5)%xy= TmpCvt3D(i  ,j+1,k-1)
                    Ppt%Cvt(6)%xy= TmpCvt3D(i-1,j-1,k  )
                    Ppt%Cvt(7)%xy= TmpCvt3D(i  ,j-1,k  )
                    Ppt%Cvt(8)%xy= TmpCvt3D(i+1,j-1,k  )
                    Ppt%Cvt(9)%xy= TmpCvt3D(i-1,j  ,k  )
                    Ppt%Cvt(10)%xy=TmpCvt3D(i  ,j  ,k  )
                    Ppt%Cvt(11)%xy=TmpCvt3D(i+1,j  ,k  )
                    Ppt%Cvt(12)%xy=TmpCvt3D(i-1,j+1,k  )
                    Ppt%Cvt(13)%xy=TmpCvt3D(i  ,j+1,k  )
                    Ppt%Cvt(14)%xy=TmpCvt3D(i+1,j+1,k  )
                    Ppt%Cvt(15)%xy=TmpCvt3D(i  ,j-1,k+1)
                    Ppt%Cvt(16)%xy=TmpCvt3D(i-1,j  ,k+1)
                    Ppt%Cvt(17)%xy=TmpCvt3D(i  ,j  ,k+1)
                    Ppt%Cvt(18)%xy=TmpCvt3D(i+1,j  ,k+1)
                    Ppt%Cvt(19)%xy=TmpCvt3D(i  ,j+1,k+1)
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%xy,TmpCvt3D,ierr)

                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%xz,TmpCvt3D,ierr)
                    Ppt%Cvt(1)%xz= TmpCvt3D(i  ,j-1,k-1)
                    Ppt%Cvt(2)%xz= TmpCvt3D(i-1,j  ,k-1)
                    Ppt%Cvt(3)%xz= TmpCvt3D(i  ,j  ,k-1)
                    Ppt%Cvt(4)%xz= TmpCvt3D(i+1,j  ,k-1)
                    Ppt%Cvt(5)%xz= TmpCvt3D(i  ,j+1,k-1)
                    Ppt%Cvt(6)%xz= TmpCvt3D(i-1,j-1,k  )
                    Ppt%Cvt(7)%xz= TmpCvt3D(i  ,j-1,k  )
                    Ppt%Cvt(8)%xz= TmpCvt3D(i+1,j-1,k  )
                    Ppt%Cvt(9)%xz= TmpCvt3D(i-1,j  ,k  )
                    Ppt%Cvt(10)%xz=TmpCvt3D(i  ,j  ,k  )
                    Ppt%Cvt(11)%xz=TmpCvt3D(i+1,j  ,k  )
                    Ppt%Cvt(12)%xz=TmpCvt3D(i-1,j+1,k  )
                    Ppt%Cvt(13)%xz=TmpCvt3D(i  ,j+1,k  )
                    Ppt%Cvt(14)%xz=TmpCvt3D(i+1,j+1,k  )
                    Ppt%Cvt(15)%xz=TmpCvt3D(i  ,j-1,k+1)
                    Ppt%Cvt(16)%xz=TmpCvt3D(i-1,j  ,k+1)
                    Ppt%Cvt(17)%xz=TmpCvt3D(i  ,j  ,k+1)
                    Ppt%Cvt(18)%xz=TmpCvt3D(i+1,j  ,k+1)
                    Ppt%Cvt(19)%xz=TmpCvt3D(i  ,j+1,k+1)
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%xz,TmpCvt3D,ierr)

                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%yz,TmpCvt3D,ierr)
                    Ppt%Cvt(1)%yz= TmpCvt3D(i  ,j-1,k-1)
                    Ppt%Cvt(2)%yz= TmpCvt3D(i-1,j  ,k-1)
                    Ppt%Cvt(3)%yz= TmpCvt3D(i  ,j  ,k-1)
                    Ppt%Cvt(4)%yz= TmpCvt3D(i+1,j  ,k-1)
                    Ppt%Cvt(5)%yz= TmpCvt3D(i  ,j+1,k-1)
                    Ppt%Cvt(6)%yz= TmpCvt3D(i-1,j-1,k  )
                    Ppt%Cvt(7)%yz= TmpCvt3D(i  ,j-1,k  )
                    Ppt%Cvt(8)%yz= TmpCvt3D(i+1,j-1,k  )
                    Ppt%Cvt(9)%yz= TmpCvt3D(i-1,j  ,k  )
                    Ppt%Cvt(10)%yz=TmpCvt3D(i  ,j  ,k  )
                    Ppt%Cvt(11)%yz=TmpCvt3D(i+1,j  ,k  )
                    Ppt%Cvt(12)%yz=TmpCvt3D(i-1,j+1,k  )
                    Ppt%Cvt(13)%yz=TmpCvt3D(i  ,j+1,k  )
                    Ppt%Cvt(14)%yz=TmpCvt3D(i+1,j+1,k  )
                    Ppt%Cvt(15)%yz=TmpCvt3D(i  ,j-1,k+1)
                    Ppt%Cvt(16)%yz=TmpCvt3D(i-1,j  ,k+1)
                    Ppt%Cvt(17)%yz=TmpCvt3D(i  ,j  ,k+1)
                    Ppt%Cvt(18)%yz=TmpCvt3D(i+1,j  ,k+1)
                    Ppt%Cvt(19)%yz=TmpCvt3D(i  ,j+1,k+1)
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%yz,TmpCvt3D,ierr)

                    CALL TargetFullTensor(Ppt%Cvt(1))
                    CALL TargetFullTensor(Ppt%Cvt(2))
                    CALL TargetFullTensor(Ppt%Cvt(3))
                    CALL TargetFullTensor(Ppt%Cvt(4))
                    CALL TargetFullTensor(Ppt%Cvt(5))
                    CALL TargetFullTensor(Ppt%Cvt(6))
                    CALL TargetFullTensor(Ppt%Cvt(7))
                    CALL TargetFullTensor(Ppt%Cvt(8))
                    CALL TargetFullTensor(Ppt%Cvt(9))
                    CALL TargetFullTensor(Ppt%Cvt(10))
                    CALL TargetFullTensor(Ppt%Cvt(11))
                    CALL TargetFullTensor(Ppt%Cvt(12))
                    CALL TargetFullTensor(Ppt%Cvt(13))
                    CALL TargetFullTensor(Ppt%Cvt(14))
                    CALL TargetFullTensor(Ppt%Cvt(15))
                    CALL TargetFullTensor(Ppt%Cvt(16))
                    CALL TargetFullTensor(Ppt%Cvt(17))
                    CALL TargetFullTensor(Ppt%Cvt(18))
                    CALL TargetFullTensor(Ppt%Cvt(19))

                ELSE ! 2D mesh
                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%xx,TmpCvt2D,ierr)
                    Ppt%Cvt(6)%xx= TmpCvt2D(i-1,j-1)
                    Ppt%Cvt(7)%xx= TmpCvt2D(i  ,j-1)
                    Ppt%Cvt(8)%xx= TmpCvt2D(i+1,j-1)
                    Ppt%Cvt(9)%xx= TmpCvt2D(i-1,j  )
                    Ppt%Cvt(10)%xx=TmpCvt2D(i  ,j  )
                    Ppt%Cvt(11)%xx=TmpCvt2D(i+1,j  )
                    Ppt%Cvt(12)%xx=TmpCvt2D(i-1,j+1)
                    Ppt%Cvt(13)%xx=TmpCvt2D(i  ,j+1)
                    Ppt%Cvt(14)%xx=TmpCvt2D(i+1,j+1)
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%xx,TmpCvt2D,ierr)

                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%yy,TmpCvt2D,ierr)
                    Ppt%Cvt(6)%yy= TmpCvt2D(i-1,j-1)
                    Ppt%Cvt(7)%yy= TmpCvt2D(i  ,j-1)
                    Ppt%Cvt(8)%yy= TmpCvt2D(i+1,j-1)
                    Ppt%Cvt(9)%yy= TmpCvt2D(i-1,j  )
                    Ppt%Cvt(10)%yy=TmpCvt2D(i  ,j  )
                    Ppt%Cvt(11)%yy=TmpCvt2D(i+1,j  )
                    Ppt%Cvt(12)%yy=TmpCvt2D(i-1,j+1)
                    Ppt%Cvt(13)%yy=TmpCvt2D(i  ,j+1)
                    Ppt%Cvt(14)%yy=TmpCvt2D(i+1,j+1)
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%yy,TmpCvt2D,ierr)

                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%zz,TmpCvt2D,ierr)
                    Ppt%Cvt(6)%zz= TmpCvt2D(i-1,j-1)
                    Ppt%Cvt(7)%zz= TmpCvt2D(i  ,j-1)
                    Ppt%Cvt(8)%zz= TmpCvt2D(i+1,j-1)
                    Ppt%Cvt(9)%zz= TmpCvt2D(i-1,j  )
                    Ppt%Cvt(10)%zz=TmpCvt2D(i  ,j  )
                    Ppt%Cvt(11)%zz=TmpCvt2D(i+1,j  )
                    Ppt%Cvt(12)%zz=TmpCvt2D(i-1,j+1)
                    Ppt%Cvt(13)%zz=TmpCvt2D(i  ,j+1)
                    Ppt%Cvt(14)%zz=TmpCvt2D(i+1,j+1)
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%zz,TmpCvt2D,ierr)

                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%xy,TmpCvt2D,ierr)
                    Ppt%Cvt(6)%xy= TmpCvt2D(i-1,j-1)
                    Ppt%Cvt(7)%xy= TmpCvt2D(i  ,j-1)
                    Ppt%Cvt(8)%xy= TmpCvt2D(i+1,j-1)
                    Ppt%Cvt(9)%xy= TmpCvt2D(i-1,j  )
                    Ppt%Cvt(10)%xy=TmpCvt2D(i  ,j  )
                    Ppt%Cvt(11)%xy=TmpCvt2D(i+1,j  )
                    Ppt%Cvt(12)%xy=TmpCvt2D(i-1,j+1)
                    Ppt%Cvt(13)%xy=TmpCvt2D(i  ,j+1)
                    Ppt%Cvt(14)%xy=TmpCvt2D(i+1,j+1)
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%xy,TmpCvt2D,ierr)

                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%xz,TmpCvt2D,ierr)
                    Ppt%Cvt(6)%xz= TmpCvt2D(i-1,j-1)
                    Ppt%Cvt(7)%xz= TmpCvt2D(i  ,j-1)
                    Ppt%Cvt(8)%xz= TmpCvt2D(i+1,j-1)
                    Ppt%Cvt(9)%xz= TmpCvt2D(i-1,j  )
                    Ppt%Cvt(10)%xz=TmpCvt2D(i  ,j  )
                    Ppt%Cvt(11)%xz=TmpCvt2D(i+1,j  )
                    Ppt%Cvt(12)%xz=TmpCvt2D(i-1,j+1)
                    Ppt%Cvt(13)%xz=TmpCvt2D(i  ,j+1)
                    Ppt%Cvt(14)%xz=TmpCvt2D(i+1,j+1)
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%xz,TmpCvt2D,ierr)

                    CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%yz,TmpCvt2D,ierr)
                    Ppt%Cvt(6)%yz= TmpCvt2D(i-1,j-1)
                    Ppt%Cvt(7)%yz= TmpCvt2D(i  ,j-1)
                    Ppt%Cvt(8)%yz= TmpCvt2D(i+1,j-1)
                    Ppt%Cvt(9)%yz= TmpCvt2D(i-1,j  )
                    Ppt%Cvt(10)%yz=TmpCvt2D(i  ,j  )
                    Ppt%Cvt(11)%yz=TmpCvt2D(i+1,j  )
                    Ppt%Cvt(12)%yz=TmpCvt2D(i-1,j+1)
                    Ppt%Cvt(13)%yz=TmpCvt2D(i  ,j+1)
                    Ppt%Cvt(14)%yz=TmpCvt2D(i+1,j+1)
                    CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%yz,TmpCvt2D,ierr)

                    CALL TargetFullTensor(Ppt%Cvt(6))
                    CALL TargetFullTensor(Ppt%Cvt(7))
                    CALL TargetFullTensor(Ppt%Cvt(8))
                    CALL TargetFullTensor(Ppt%Cvt(9))
                    CALL TargetFullTensor(Ppt%Cvt(10))
                    CALL TargetFullTensor(Ppt%Cvt(11))
                    CALL TargetFullTensor(Ppt%Cvt(12))
                    CALL TargetFullTensor(Ppt%Cvt(13))
                    CALL TargetFullTensor(Ppt%Cvt(14))

                END IF      
            ELSE
                PRINT*,"ERROR: The properties has to be defined by cell or by zone"
            END IF  
        END IF
        IF (k.EQ.0) Ppt%Cvt(1:5)=TensorZero
        IF (j.EQ.0) THEN
            Ppt%Cvt(1)=TensorZero;Ppt%Cvt(6:8)=TensorZero;Ppt%Cvt(15)=TensorZero
        END IF
        IF (i.EQ.0) THEN
            Ppt%Cvt(2)=TensorZero;Ppt%Cvt(6)=TensorZero;Ppt%Cvt(9)=TensorZero;Ppt%Cvt(12)=TensorZero;Ppt%Cvt(16)=TensorZero
        END IF
        IF (i.EQ.(widthG(1)-1)) THEN
            Ppt%Cvt(4)=TensorZero;Ppt%Cvt(8)=TensorZero;Ppt%Cvt(11)=TensorZero;Ppt%Cvt(14)=TensorZero;Ppt%Cvt(18)=TensorZero
        END IF
        IF (j.EQ.(widthG(2)-1)) THEN
            Ppt%Cvt(5)=TensorZero;Ppt%Cvt(12:14)=TensorZero;Ppt%Cvt(19)=TensorZero
        END IF
        IF (k.EQ.(widthG(3)-1)) Ppt%Cvt(15:19)=TensorZero
    END IF

END SUBROUTINE GetLocalConductivity

SUBROUTINE DestroyProperties(PptFld,ierr)

    USE ANISOFLOW_Types, ONLY : PropertiesField,InputTypeVar
    USE ANISOFLOW_Interface, ONLY : GetVerbose

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    TYPE(PropertiesField),INTENT(INOUT) :: PptFld

    PetscBool                           :: Verbose
    CHARACTER(LEN=200)                  :: EventName,ClassName
    PetscLogEvent                       :: Event
    PetscClassId                        :: ClassID
    PetscLogDouble                      :: EventFlops=0.d0

    ClassName="Property"
    CALL PetscClassIdRegister(ClassName,ClassID,ierr)
    EventName="DestroyProperties"
    CALL PetscLogEventRegister(EventName,ClassID,Event,ierr)
    CALL PetscLogEventBegin(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Inizialited\n",ierr)

!     IF (PptFld%Cvt%DefinedByZones) THEN
!         IF (ALLOCATED(PptFld%Cvt%CvtZone)) DEALLOCATE(PptFld%Cvt%CvtZone)
!         CALL VecDestroy(PptFld%PptType,ierr)
!     ELSE
!         CALL VecDestroy(PptFld%Cvt%Cvt,ierr)
!     END IF

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Finalized\n",ierr)
    
    CALL PetscLogFlops(EventFlops,ierr)
    CALL PetscLogEventEnd(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

END SUBROUTINE DestroyProperties

END MODULE ANISOFLOW_Properties