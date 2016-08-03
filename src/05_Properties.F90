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
    PetscLogDouble                      :: EventFlops=0.d0
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
    USE ANISOFLOW_Interface, ONLY : GetVerbose,GetInputDir,GetInputFilePptByZone
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
    Vec,INTENT(OUT)                     :: PptZoneID_Local
    PetscInt,INTENT(OUT)                :: DefinedByPptZones

    PetscInt                            :: widthG(3),u,ValI,i
    PetscMPIInt                         :: process
    PetscReal                           :: ValR
    CHARACTER(LEN=200)                  :: InputFilePptByZone,Route,ViewName,InputDir!,EventName
    PetscBool                           :: InputFilePptByZoneFlg,Verbose
    Vec                                 :: PptZoneID_Global

    PARAMETER(u=01)

    CALL GetVerbose(Verbose,ierr)

    CALL GetInputDir(InputDir,ierr)
    CALL GetInputFilePptByZone(InputFilePptByZone,InputFilePptByZoneFlg,ierr)

    IF (InputFilePptByZoneFlg) THEN

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
                            & "[ERROR] Input_file_ppt_by_zones entry only can contain"  &
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

        ViewName="ANISOFLOW_PptZoneID"
!         EventName="GetConductivity"
        ! Crear un visualizador para esto
!         CALL ViewConductivity(PptZoneID,ViewName,EventName,ierr)

        CALL DMGlobalToLocalBegin(Gmtry%DataMngr,PptZoneID_Global,INSERT_VALUES,PptZoneID_Local,ierr)
        CALL DMGlobalToLocalEnd(Gmtry%DataMngr,PptZoneID_Global,INSERT_VALUES,PptZoneID_Local,ierr)
        CALL VecDestroy(PptZoneID_Global,ierr)
        DefinedByPptZones=1
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
    ELSE
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "[ERROR] Conductivity InputType wrong\n",ierr)
        STOP
    END IF

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Finalized\n",ierr)
    
    CALL PetscLogFlops(EventFlops,ierr)
    CALL PetscLogEventEnd(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

END SUBROUTINE GetConductivity

SUBROUTINE GetCvtZoneID(Gmtry,PptFld,CvtZoneID_Local,DefinedByPptZones,DefinedByCvtZones,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,PropertiesField,TargPetscVec
    USE ANISOFLOW_Interface, ONLY : GetVerbose,GetInputDir,GetInputFilePptByZone,GetInputFileCvtByZone
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
    PetscInt,INTENT(OUT)                :: DefinedByPptZones
    PetscInt,INTENT(OUT)                :: DefinedByCvtZones

    PetscInt                            :: widthG(3),u,ValI,i
    PetscMPIInt                         :: process
    PetscReal                           :: ValR
    CHARACTER(LEN=200)                  :: InputFilePptByZone,InputFileCvtByZone,Route,ViewName,InputDir!,EventName
    PetscBool                           :: InputFilePptByZoneFlg,InputFileCvtByZoneFlg,Verbose
    Vec                                 :: CvtZoneID_Global
    Vec,POINTER                         :: ZoneID_tmp

    PARAMETER(u=01)

    CALL GetVerbose(Verbose,ierr)

    CALL GetInputDir(InputDir,ierr)
    CALL GetInputFilePptByZone(InputFilePptByZone,InputFilePptByZoneFlg,ierr)
    CALL GetInputFileCvtByZone(InputFileCvtByZone,InputFileCvtByZoneFlg,ierr)

    IF (InputFileCvtByZoneFlg) THEN

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

        ViewName="ANISOFLOW_CvtZoneID"
!         EventName="GetConductivity"
        ! Crear un visualizador para esto
!         CALL ViewConductivity(CvtZoneID,ViewName,EventName,ierr)

        CALL DMGlobalToLocalBegin(Gmtry%DataMngr,CvtZoneID_Global,INSERT_VALUES,CvtZoneID_Local,ierr)
        CALL DMGlobalToLocalEnd(Gmtry%DataMngr,CvtZoneID_Global,INSERT_VALUES,CvtZoneID_Local,ierr)
        CALL VecDestroy(CvtZoneID_Global,ierr)

        DefinedByCvtZones=1
    ELSEIF (InputFilePptByZoneFlg) THEN
        ZoneID_tmp => TargPetscVec(PptFld%ZoneID)
        CvtZoneID_Local = ZoneID_tmp
        DefinedByPptZones=1
    ELSE
        PRINT*,"ERROR en CvtZone"! TODO: arreglar este mensage
        STOP
    END IF

END SUBROUTINE GetCvtZoneID

SUBROUTINE GetConductivity_1(Gmtry,PptFld,Cvt,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,ConductivityField,PropertiesField,TargetFullTensor
    USE ANISOFLOW_Interface, ONLY : GetVerbose,GetInputDir,GetInputFileCvt
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
    TYPE(ConductivityField),INTENT(OUT) :: Cvt

    PetscInt                            :: u,i,j,CvtLen
    PetscReal                           :: ValR
    PetscMPIInt                         :: process
    CHARACTER(LEN=200)                  :: InputDir,InputFileCvt!,ViewName!,EventName
    CHARACTER(LEN=200)                  :: Route
    CHARACTER(LEN=13)                   :: CvtKind
    PetscBool                           :: Verbose

    PARAMETER(u=01)

    CALL GetVerbose(Verbose,ierr)

    CALL GetCvtZoneID(Gmtry,PptFld,Cvt%ZoneID,Cvt%DefinedByPptZones,Cvt%DefinedByCvtZones,ierr)
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
            READ(u,'(A13)')CvtKind
            IF (CvtKind.EQ."k anisotropic") THEN 
                READ(u,*)Cvt%Zone(i)%xx, &
                    & Cvt%Zone(i)%yy,Cvt%Zone(i)%zz
                Cvt%Zone(i)%xy=0.0
                Cvt%Zone(i)%xz=0.0
                Cvt%Zone(i)%yz=0.0
            ELSE IF (CvtKind.EQ."k isotropic  ") THEN
                READ(u,*)Cvt%Zone(i)%xx
                Cvt%Zone(i)%yy=Cvt%Zone(i)%xx
                Cvt%Zone(i)%zz=Cvt%Zone(i)%xx
                Cvt%Zone(i)%xy=0.0
                Cvt%Zone(i)%xz=0.0
                Cvt%Zone(i)%yz=0.0
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
    PetscReal                           :: ValR
    PetscMPIInt                         :: process
    CHARACTER(LEN=200)                  :: InputDir,InputFileCvt
    CHARACTER(LEN=200)                  :: Route!,ViewName,EventName
    PetscBool                           :: Verbose
    Vec                                 :: Cell_Global

    PARAMETER(u=01)

    CALL GetVerbose(Verbose,ierr)

    CALL GetInputDir(InputDir,ierr)
    CALL GetInputFileCvt(InputFileCvt,ierr)

    Cvt%DefinedByCell=1
!     IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetProrperties Event] Conductivity Field is stored by Block\n",ierr)

    CALL DMCreateLocalVector(Gmtry%DataMngr,Cvt%Cell,ierr)
    CALL DMCreateGlobalVector(Gmtry%DataMngr,Cell_Global,ierr)

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
            CALL VecSetValue(Cell_Global,i-1,ValR,INSERT_VALUES,ierr)
        END DO
        CLOSE(u)
    END IF

    CALL VecAssemblyBegin(Cell_Global,ierr)
    CALL VecAssemblyEnd(Cell_Global,ierr)

!     ViewName="ANISOFLOW_Cvt"
!     EventName="GetConductivity"
!     CALL ViewConductivity(Cvt%Cell,ViewName,EventName,ierr)

    CALL DMGlobalToLocalBegin(Gmtry%DataMngr,Cell_Global,INSERT_VALUES,Cvt%Cell,ierr)
    CALL DMGlobalToLocalEnd(Gmtry%DataMngr,Cell_Global,INSERT_VALUES,Cvt%Cell,ierr)
    CALL VecDestroy(Cell_Global,ierr)

!     IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetProrperties Stage] Conductivity Field was satisfactorily created\n",ierr)

END SUBROUTINE GetConductivity_2


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

SUBROUTINE GetStoZoneID(Gmtry,PptFld,StoZoneID_Local,DefinedByPptZones,DefinedByStoZones,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,PropertiesField,TargPetscVec
    USE ANISOFLOW_Interface, ONLY : GetVerbose,GetInputDir,GetInputFilePptByZone,GetInputFileStoByZone
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
    Vec,INTENT(OUT)                     :: StoZoneID_Local
    PetscInt,INTENT(OUT)                :: DefinedByPptZones,DefinedByStoZones

    PetscInt                            :: widthG(3),u,ValI,i
    PetscMPIInt                         :: process
    PetscReal                           :: ValR
    CHARACTER(LEN=200)                  :: InputFilePptByZone,InputFileStoByZone,Route,ViewName,InputDir!,EventName
    PetscBool                           :: InputFilePptByZoneFlg,InputFileStoByZoneFlg,Verbose
    Vec,POINTER                         :: ZoneID_tmp
    Vec                                 :: StoZoneID_Global

    PARAMETER(u=01)

    CALL GetVerbose(Verbose,ierr)

    CALL GetInputDir(InputDir,ierr)
    CALL GetInputFilePptByZone(InputFilePptByZone,InputFilePptByZoneFlg,ierr)
    CALL GetInputFileStoByZone(InputFileStoByZone,InputFileStoByZoneFlg,ierr)


    IF (InputFileStoByZoneFlg) THEN

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

        ViewName="ANISOFLOW_StoZoneID"
!         EventName="GetConductivity"
        ! Crear un visualizador para esto
!         CALL ViewConductivity(StoZoneID,ViewName,EventName,ierr)

        CALL DMGlobalToLocalBegin(Gmtry%DataMngr,StoZoneID_Global,INSERT_VALUES,StoZoneID_Local,ierr)
        CALL DMGlobalToLocalEnd(Gmtry%DataMngr,StoZoneID_Global,INSERT_VALUES,StoZoneID_Local,ierr)
        CALL VecDestroy(StoZoneID_Global,ierr)

        DefinedByStoZones=1
    ELSEIF (InputFilePptByZoneFlg) THEN
        ZoneID_tmp => TargPetscVec(PptFld%ZoneID)
        StoZoneID_Local = ZoneID_tmp
        DefinedByPptZones=1
    ELSE
        PRINT*,"ERROR en StoZone"! TODO: arreglar este mensage
        STOP
    END IF

END SUBROUTINE GetStoZoneID

SUBROUTINE GetStorage_1(Gmtry,PptFld,Sto,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,SpecificStorageField,PropertiesField,TargetFullTensor
    USE ANISOFLOW_Interface, ONLY : GetVerbose,GetInputDir,GetInputFileSto
!     USE ANISOFLOW_View, ONLY : ViewStorage
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

    CALL GetStoZoneID(Gmtry,PptFld,Sto%ZoneID,Sto%DefinedByPptZones,Sto%DefinedByStoZones,ierr)
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
    USE ANISOFLOW_Interface, !ONLY : GetVerbose,GetInputDir,GetInputFileSto

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

    PetscInt                :: widthL(3),widthG(3),corn(3),ValI,i,j,k
    PetscReal,POINTER       :: TmpStoZone(:),TmpStoZoneID3D(:,:,:),TmpStoZoneID2D(:,:),TmpStoCell3D(:,:,:),TmpStoCell2D(:,:)

    CALL DMDAGetInfo(Gmtry%DataMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
        & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,ierr)

    IF ((Sto%DefinedByPptZones.EQ.0).AND.(Sto%DefinedByStoZones.EQ.0)) THEN
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
                    TmpStoZoneID3D(i,j,k)=TmpStoZone(ValI)
                ELSE
                    ValI=INT(TmpStoZoneID2D(i,j))
                    TmpStoZoneID2D(i,j)=TmpStoZone(ValI)
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

    IF (Sto%DefinedByStoZones.EQ.1) THEN
        CALL VecDestroy(Sto%ZoneID,ierr)
        Sto%DefinedByStoZones=0
    ELSE IF (Sto%DefinedByPptZones.EQ.1) THEN
        !NULLIFY(Sto%ZoneID) ! Because I used a very tricky way to use it as a pointer, I can't nullify it as a usual pointer but it's supposed it doesn't will get in troubles the code.
        Sto%DefinedByPptZones=0
    END IF

    Sto%DefinedByCell=1

!     CALL VecView(Sto%Cell,PETSC_VIEWER_STDOUT_WORLD,ierr)

END SUBROUTINE StorageZoneToCell

SUBROUTINE GetStorage_2(Gmtry,PptFld,Sto,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,PropertiesField,SpecificStorageField
    USE ANISOFLOW_Interface, ONLY : GetVerbose,GetInputDir,GetInputFileSto
!     USE ANISOFLOW_View, ONLY : ViewStorage
    
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
    CHARACTER(LEN=200)                  :: Route!,ViewName,EventName
    PetscBool                           :: Verbose

    PARAMETER(u=01)

    CALL GetVerbose(Verbose,ierr)

    CALL GetInputDir(InputDir,ierr)
    CALL GetInputFileSto(InputFileSto,ierr)

    Sto%DefinedByCell=1
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

!     ViewName="ANISOFLOW_Sto"
!     EventName="GetStorage"
!     CALL ViewStorage(Sto%Cell,ViewName,EventName,ierr)

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

    USE ANISOFLOW_Types, ONLY : Geometry,PropertiesField,Property,TargetFullTensor,InputTypeVar
    USE ANISOFLOW_Interface, ONLY : GetVerbose,GetInputType

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

    TYPE(InputTypeVar)                  :: InputType
    PetscReal,POINTER                   :: TmpCvt3D(:,:,:),TmpCvt2D(:,:)
    PetscReal                           :: ValR(2)
    PetscInt                            :: ValI(2),widthG(3),i,j,k

    CALL GetInputType(InputType,ierr)

    i=Ppt%Pstn%i
    j=Ppt%Pstn%j
    k=Ppt%Pstn%k

    ! It gets the global size from the geometry data manager.
    CALL DMDAGetInfo(Gmtry%DataMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
        & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,ierr)

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

    IF ((PptFld%Cvt%DefinedByCvtZones.EQ.1).OR.(PptFld%Cvt%DefinedByPptZones.EQ.1)) THEN
        IF ((Ppt%StnclTplgy(ValI(1)).EQ.1).OR.(Ppt%StnclTplgy(ValI(1)).EQ.3).OR.(Ppt%StnclTplgy(ValI(1)).EQ.4).OR.(Ppt%StnclTplgy(ValI(1)).EQ.5)) THEN ! Only get properties for active blocks
            IF (widthG(3).NE.1) THEN
                CALL DMDAVecGetArrayreadF90(Gmtry%DataMngr,PptFld%Cvt%ZoneID,TmpCvt3D,ierr)
                ValI(1)=INT(TmpCvt3D(i,j,k))
            ELSE
                CALL DMDAVecGetArrayreadF90(Gmtry%DataMngr,PptFld%Cvt%ZoneID,TmpCvt2D,ierr)
                ValI(1)=INT(TmpCvt2D(i,j))
            END IF                

            ! Backward on x
            IF (widthG(3).NE.1) THEN
                ValI(2)=INT(TmpCvt3D(i-1,j,k))
            ELSE
                ValI(2)=INT(TmpCvt2D(i-1,j))
            END IF
            Ppt%CvtBx%xx=Armonic(PptFld%Cvt%Zone(ValI(1))%xx,PptFld%Cvt%Zone(ValI(2))%xx)
            Ppt%CvtBx%yy=Armonic(PptFld%Cvt%Zone(ValI(1))%yy,PptFld%Cvt%Zone(ValI(2))%yy)
            Ppt%CvtBx%zz=Armonic(PptFld%Cvt%Zone(ValI(1))%zz,PptFld%Cvt%Zone(ValI(2))%zz)
            Ppt%CvtBx%xy=Armonic(PptFld%Cvt%Zone(ValI(1))%xy,PptFld%Cvt%Zone(ValI(2))%xy)
            Ppt%CvtBx%xz=Armonic(PptFld%Cvt%Zone(ValI(1))%xz,PptFld%Cvt%Zone(ValI(2))%xz)
            Ppt%CvtBx%yz=Armonic(PptFld%Cvt%Zone(ValI(1))%yz,PptFld%Cvt%Zone(ValI(2))%yz)
            CALL TargetFullTensor(Ppt%CvtBx)

            ! Forward on x
            IF (widthG(3).NE.1) THEN
                ValI(2)=INT(TmpCvt3D(i+1,j,k))
            ELSE
                ValI(2)=INT(TmpCvt2D(i+1,j))
            END IF
            Ppt%CvtFx%xx=Armonic(PptFld%Cvt%Zone(ValI(1))%xx,PptFld%Cvt%Zone(ValI(2))%xx)
            Ppt%CvtFx%yy=Armonic(PptFld%Cvt%Zone(ValI(1))%yy,PptFld%Cvt%Zone(ValI(2))%yy)
            Ppt%CvtFx%zz=Armonic(PptFld%Cvt%Zone(ValI(1))%zz,PptFld%Cvt%Zone(ValI(2))%zz)
            Ppt%CvtFx%xy=Armonic(PptFld%Cvt%Zone(ValI(1))%xy,PptFld%Cvt%Zone(ValI(2))%xy)
            Ppt%CvtFx%xz=Armonic(PptFld%Cvt%Zone(ValI(1))%xz,PptFld%Cvt%Zone(ValI(2))%xz)
            Ppt%CvtFx%yz=Armonic(PptFld%Cvt%Zone(ValI(1))%yz,PptFld%Cvt%Zone(ValI(2))%yz)
            CALL TargetFullTensor(Ppt%CvtFx)

            ! Backward on y
            IF (widthG(3).NE.1) THEN
                ValI(2)=INT(TmpCvt3D(i,j-1,k))
            ELSE
                ValI(2)=INT(TmpCvt2D(i-1,j))
            END IF
            Ppt%CvtBy%xx=Armonic(PptFld%Cvt%Zone(ValI(1))%xx,PptFld%Cvt%Zone(ValI(2))%xx)
            Ppt%CvtBy%yy=Armonic(PptFld%Cvt%Zone(ValI(1))%yy,PptFld%Cvt%Zone(ValI(2))%yy)
            Ppt%CvtBy%zz=Armonic(PptFld%Cvt%Zone(ValI(1))%zz,PptFld%Cvt%Zone(ValI(2))%zz)
            Ppt%CvtBy%xy=Armonic(PptFld%Cvt%Zone(ValI(1))%xy,PptFld%Cvt%Zone(ValI(2))%xy)
            Ppt%CvtBy%xz=Armonic(PptFld%Cvt%Zone(ValI(1))%xz,PptFld%Cvt%Zone(ValI(2))%xz)
            Ppt%CvtBy%yz=Armonic(PptFld%Cvt%Zone(ValI(1))%yz,PptFld%Cvt%Zone(ValI(2))%yz)
            CALL TargetFullTensor(Ppt%CvtBy)

            ! Forward on y
            IF (widthG(3).NE.1) THEN
                ValI(2)=INT(TmpCvt3D(i,j+1,k))
            ELSE
                ValI(2)=INT(TmpCvt2D(i,j+1))
            END IF
            Ppt%CvtFy%xx=Armonic(PptFld%Cvt%Zone(ValI(1))%xx,PptFld%Cvt%Zone(ValI(2))%xx)
            Ppt%CvtFy%yy=Armonic(PptFld%Cvt%Zone(ValI(1))%yy,PptFld%Cvt%Zone(ValI(2))%yy)
            Ppt%CvtFy%zz=Armonic(PptFld%Cvt%Zone(ValI(1))%zz,PptFld%Cvt%Zone(ValI(2))%zz)
            Ppt%CvtFy%xy=Armonic(PptFld%Cvt%Zone(ValI(1))%xy,PptFld%Cvt%Zone(ValI(2))%xy)
            Ppt%CvtFy%xz=Armonic(PptFld%Cvt%Zone(ValI(1))%xz,PptFld%Cvt%Zone(ValI(2))%xz)
            Ppt%CvtFy%yz=Armonic(PptFld%Cvt%Zone(ValI(1))%yz,PptFld%Cvt%Zone(ValI(2))%yz)
            CALL TargetFullTensor(Ppt%CvtFy)

            ! Backward on z
            IF (widthG(3).NE.1) THEN
                ValI(2)=INT(TmpCvt3D(i,j,k-1))
                Ppt%CvtBz%xx=Armonic(PptFld%Cvt%Zone(ValI(1))%xx,PptFld%Cvt%Zone(ValI(2))%xx)
                Ppt%CvtBz%yy=Armonic(PptFld%Cvt%Zone(ValI(1))%yy,PptFld%Cvt%Zone(ValI(2))%yy)
                Ppt%CvtBz%zz=Armonic(PptFld%Cvt%Zone(ValI(1))%zz,PptFld%Cvt%Zone(ValI(2))%zz)
                Ppt%CvtBz%xy=Armonic(PptFld%Cvt%Zone(ValI(1))%xy,PptFld%Cvt%Zone(ValI(2))%xy)
                Ppt%CvtBz%xz=Armonic(PptFld%Cvt%Zone(ValI(1))%xz,PptFld%Cvt%Zone(ValI(2))%xz)
                Ppt%CvtBz%yz=Armonic(PptFld%Cvt%Zone(ValI(1))%yz,PptFld%Cvt%Zone(ValI(2))%yz)
                CALL TargetFullTensor(Ppt%CvtBz)
            END IF


            ! Forward on y
            IF (widthG(3).NE.1) THEN
                ValI(2)=INT(TmpCvt3D(i,j,k+1))
                Ppt%CvtFz%xx=Armonic(PptFld%Cvt%Zone(ValI(1))%xx,PptFld%Cvt%Zone(ValI(2))%xx)
                Ppt%CvtFz%yy=Armonic(PptFld%Cvt%Zone(ValI(1))%yy,PptFld%Cvt%Zone(ValI(2))%yy)
                Ppt%CvtFz%zz=Armonic(PptFld%Cvt%Zone(ValI(1))%zz,PptFld%Cvt%Zone(ValI(2))%zz)
                Ppt%CvtFz%xy=Armonic(PptFld%Cvt%Zone(ValI(1))%xy,PptFld%Cvt%Zone(ValI(2))%xy)
                Ppt%CvtFz%xz=Armonic(PptFld%Cvt%Zone(ValI(1))%xz,PptFld%Cvt%Zone(ValI(2))%xz)
                Ppt%CvtFz%yz=Armonic(PptFld%Cvt%Zone(ValI(1))%yz,PptFld%Cvt%Zone(ValI(2))%yz)
                CALL TargetFullTensor(Ppt%CvtFz)
            END IF

            IF (widthG(3).NE.1) THEN
                CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%ZoneID,TmpCvt3D,ierr)
            ELSE
                CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%ZoneID,TmpCvt2D,ierr)
            END IF
        END IF
    ELSE IF (Pptfld%Cvt%DefinedByCell.EQ.1) THEN
        IF ((Ppt%StnclTplgy(ValI(1)).EQ.1).OR.(Ppt%StnclTplgy(ValI(1)).EQ.3).OR.(Ppt%StnclTplgy(ValI(1)).EQ.4).OR.(Ppt%StnclTplgy(ValI(1)).EQ.5)) THEN

            IF (widthG(3).NE.1) THEN
                CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%Cell,TmpCvt3D,ierr)
                ValR(1)=TmpCvt3D(i,j,k)
            ELSE
                CALL DMDAVecGetArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%Cell,TmpCvt2D,ierr)
                ValR(1)=TmpCvt2D(i,j)
            END IF      

            IF (widthG(3).NE.1) THEN
                ValR(2)=TmpCvt3D(i-1,j,k)
            ELSE
                ValR(2)=TmpCvt2D(i-1,j)
            END IF  
            Ppt%CvtBx%xx=Armonic(ValR(1),ValR(2))
            Ppt%CvtBx%xy=0.0
            Ppt%CvtBx%xz=0.0
            Ppt%CvtBx%yy=0.0
            Ppt%CvtBx%yz=0.0
            Ppt%CvtBx%zz=0.0
            CALL TargetFullTensor(Ppt%CvtBx)

            IF (widthG(3).NE.1) THEN
                ValR(2)=TmpCvt3D(i+1,j,k)
            ELSE
                ValR(2)=TmpCvt2D(i+1,j)
            END IF  
            Ppt%CvtFx%xx=Armonic(ValR(1),ValR(2))
            Ppt%CvtFx%xy=0.0
            Ppt%CvtFx%xz=0.0
            Ppt%CvtFx%yy=0.0
            Ppt%CvtFx%yz=0.0
            Ppt%CvtFx%zz=0.0
            CALL TargetFullTensor(Ppt%CvtFx)

            IF (widthG(3).NE.1) THEN
                ValR(2)=TmpCvt3D(i,j-1,k)
            ELSE
                ValR(2)=TmpCvt2D(i,j-1)
            END IF  
            Ppt%CvtBy%xx=0.0
            Ppt%CvtBy%xy=0.0
            Ppt%CvtBy%xz=0.0
            Ppt%CvtBy%yy=Armonic(ValR(1),ValR(2))
            Ppt%CvtBy%yz=0.0
            Ppt%CvtBy%zz=0.0
            CALL TargetFullTensor(Ppt%CvtBy)

            IF (widthG(3).NE.1) THEN
                ValR(2)=TmpCvt3D(i,j+1,k)
            ELSE
                ValR(2)=TmpCvt2D(i,j+1)
            END IF  
            Ppt%CvtFy%xx=0.0
            Ppt%CvtFy%xy=0.0
            Ppt%CvtFy%xz=0.0
            Ppt%CvtFy%yy=Armonic(ValR(1),ValR(2))
            Ppt%CvtFy%yz=0.0
            Ppt%CvtFy%zz=0.0
            CALL TargetFullTensor(Ppt%CvtFy)

            IF (widthG(3).NE.1) THEN
                ValR(2)=TmpCvt3D(i,j,k-1)
                Ppt%CvtBz%xx=0.0
                Ppt%CvtBz%xy=0.0
                Ppt%CvtBz%xz=0.0
                Ppt%CvtBz%yy=0.0
                Ppt%CvtBz%yz=0.0
                Ppt%CvtBz%zz=Armonic(ValR(1),ValR(2))
                CALL TargetFullTensor(Ppt%CvtBz)
            END IF  

            IF (widthG(3).NE.1) THEN
                ValR(2)=TmpCvt3D(i,j,k+1)
                Ppt%CvtFz%xx=0.0
                Ppt%CvtFz%xy=0.0
                Ppt%CvtFz%xz=0.0
                Ppt%CvtFz%yy=0.0
                Ppt%CvtFz%yz=0.0
                Ppt%CvtFz%zz=Armonic(ValR(1),ValR(2))
                CALL TargetFullTensor(Ppt%CvtFz)
            END IF  
            IF (widthG(3).NE.1) THEN
                CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%Cell,TmpCvt3D,ierr)
            ELSE
                CALL DMDAVecRestoreArrayReadF90(Gmtry%DataMngr,PptFld%Cvt%Cell,TmpCvt2D,ierr)
            END IF
        ELSE
            PRINT*,"ERROR: The properties has to be defined by cell or by zone"
        END IF  
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
!         CALL VecDestroy(PptFld%Cvt%CvtCell,ierr)
!     END IF

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Finalized\n",ierr)
    
    CALL PetscLogFlops(EventFlops,ierr)
    CALL PetscLogEventEnd(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

END SUBROUTINE DestroyProperties

END MODULE ANISOFLOW_Properties