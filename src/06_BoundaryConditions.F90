MODULE ANISOFLOW_BoundaryConditions

    USE ANISOFLOW_Interface, ONLY : GetVerbose
    USE ANISOFLOW_Types, ONLY : Geometry,BoundaryConditions

    IMPLICIT NONE

CONTAINS

SUBROUTINE GetBC(Gmtry,BCFld,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Geometry),INTENT(IN)               :: Gmtry
    TYPE(BoundaryConditions),INTENT(OUT)    :: BCFld

    TYPE(InputTypeVar)                      :: InputType
    PetscLogStage                           :: Stage
    PetscBool                               :: Verbose

    CALL PetscLogStageRegister("GetBC", stage,ierr)
    CALL PetscLogStagePush(Stage,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetBC Stage] Inizialited\n",ierr)

    CALL GetInputType(InputType,ierr)

    IF (InputType%BC.EQ.1) THEN
        CALL GetBC_1(Gmtry,BCFld,ierr)
    ELSE 
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "[ERROR] Boundary Condition InputType wrong\n",ierr)
        STOP
    END IF

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetBC Stage] Finalized\n",ierr)
    CALL PetscLogStagePop(Stage,ierr)

END SUBROUTINE GetBC

SUBROUTINE GetBC_1(Gmtry,BCFld,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Geometry),INTENT(IN)               :: Gmtry
    TYPE(BoundaryConditions),INTENT(OUT)    :: BCFld

    PetscMPIInt                             :: process
    PetscInt                                :: u,i,j,CountTimeZone,CountDirich,CountNeumman,CountCauchy
    PetscInt,ALLOCATABLE                    :: IndexDirich(:),IndexNeumman(:),IndexCauchy(:)
    PetscReal                               :: ValR,DT
    CHARACTER(LEN=200)                      :: InputFileBC,InputDir,Route
    CHARACTER(LEN=20)                       :: TextTimeZones
    CHARACTER(LEN=11)                       :: Text1TimeZone
    CHARACTER(LEN=3)                        :: Text2TimeZone,Text2Dirich,Text2Neumman,Text2Cauchy
    CHARACTER(LEN=4)                        :: TextDT
    CHARACTER(LEN=27)                       :: Text1Dirich
    CHARACTER(LEN=25)                       :: Text1Neumman
    CHARACTER(LEN=24)                       :: Text1Cauchy
    PetscBool                               :: AllocationFlag
    PetscBool                               :: Verbose

    PARAMETER(u=01)

    CALL GetVerbose(Verbose,ierr)

    CALL GetInputDir(InputDir,ierr)
    CALL GetInputFileBC(InputFileBC,ierr)
    
    CALL MPI_Comm_rank(MPI_COMM_WORLD,process,ierr)

    IF (process.EQ.0) THEN
        Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFileBC))
        OPEN(u,FILE=TRIM(Route),STATUS='OLD',ACTION='READ')
        READ(u,'(A20,I8)')TextTimeZones,BCFld%SizeTimeZone
        ! Imprimir errores de lectura de texto
    END IF

    CALL MPI_Bcast(BCFld%SizeTimeZone,1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)
        
    ALLOCATE(BCFld%TimeZone(BCFld%SizeTimeZone))
    ALLOCATE(BCFld%Dirich(BCFld%SizeTimeZone))
    ALLOCATE(BCFld%Neumman(BCFld%SizeTimeZone))
    ALLOCATE(BCFld%Cauchy(BCFld%SizeTimeZone))

    DO i=1,BCFld%SizeTimeZone
        ! TIME
        IF (process.EQ.0) THEN
            READ(u,'(A9,I8,A3,I8)')Text1TimeZone,CountTimeZone,Text2TimeZone,BCFld%TimeZone(i)%SizeTime
            ! Imprimir errores de lectura de texto
        END IF

        CALL MPI_Bcast(BCFld%TimeZone(i)%SizeTime,1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)
        AllocationFlag=ALLOCATED(BCFld%TimeZone(i)%Time)
        IF (.NOT.AllocationFlag) ALLOCATE(BCFld%TimeZone(i)%Time(BCFld%TimeZone(i)%SizeTime))

        IF (process.EQ.0) THEN
            READ(u,'(A4,F15.10)')TextDT,DT
            IF ((i.EQ.1).AND.(j.EQ.1)) THEN
                BCFld%TimeZone(1)%Time(1)=0.0
            ELSEIF (j.EQ.1) THEN
                BCFld%TimeZone(i)%Time(1)=BCFld%TimeZone(i-1)%Time(BCFld%TimeZone(i-1)%SizeTime)+DT
            ELSE
                DO j=2,BCFld%TimeZone(i)%SizeTime
                    BCFld%TimeZone(i)%Time(j)=BCFld%TimeZone(i)%Time(j-1)+DT
                END DO
            END IF
            ! DO j=1,BCFld%TimeZone(i)%SizeTime
            !     READ(u,'(F15.10)')BCFld%TimeZone(i)%Time(j)
            ! END DO
        END IF
        CALL MPI_Bcast(BCFld%TimeZone(i)%Time,BCFld%TimeZone(i)%SizeTime,MPI_INT, 0, PETSC_COMM_WORLD,ierr)
        
        ! DIRICH
        IF (process.EQ.0) THEN
            READ(u,'(A27,I8,A3,I8)')Text1Dirich,CountDirich,Text2Dirich,BCFld%SizeDirich
            ! Imprimir errores de lectura de texto
        END IF
        
        CALL MPI_Bcast(BCFld%SizeDirich,1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)
        CALL VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,BCFld%SizeDirich,BCFld%Dirich(i),ierr)
        
        AllocationFlag=ALLOCATED(IndexDirich)
        IF (.NOT.AllocationFlag) ALLOCATE(IndexDirich(BCFld%SizeDirich))

        IF (process.EQ.0) THEN
            DO j=1,BCFld%SizeDirich
                READ(u, '(I12,F15.10)')IndexDirich(j),ValR
                CALL VecSetValue(BCFld%Dirich(i),j-1,ValR,INSERT_VALUES,ierr)
            END DO
        END IF

        CALL VecAssemblyBegin(BCFld%Dirich(i),ierr)
        CALL VecAssemblyEnd(BCFld%Dirich(i),ierr)
        
        ! NEUMMAN
        IF (process.EQ.0) THEN
            READ(u,'(A25,I8,A3,I8)')Text1Neumman,CountNeumman,Text2Neumman,BCFld%SizeNeumman
            ! Imprimir errores de lectura de texto
        END IF

        CALL MPI_Bcast(BCFld%SizeNeumman,1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)
        CALL VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,BCFld%SizeNeumman,BCFld%Neumman(i),ierr)

        AllocationFlag=ALLOCATED(IndexNeumman)
        IF (.NOT.AllocationFlag) ALLOCATE(IndexNeumman(BCFld%SizeNeumman))

        IF (process.EQ.0) THEN
            DO j=1,BCFld%SizeNeumman
                READ(u, '(I12,F15.10)')IndexNeumman(j),ValR
                CALL VecSetValue(BCFld%Neumman(i),j-1,ValR,INSERT_VALUES,ierr)
            END DO
        END IF

        CALL VecAssemblyBegin(BCFld%Neumman(i),ierr)
        CALL VecAssemblyEnd(BCFld%Neumman(i),ierr)

        ! CAUCHY
        IF (process.EQ.0) THEN
            READ(u,'(A24,I8,A3,I8)')Text1Cauchy,CountCauchy,Text2Cauchy,BCFld%SizeCauchy
            ! Imprimir errores de lectura de texto
        END IF

        CALL MPI_Bcast(BCFld%SizeCauchy,1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)
        CALL VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,BCFld%SizeCauchy,BCFld%Cauchy(i),ierr)

        AllocationFlag=ALLOCATED(IndexCauchy)
        IF (.NOT.AllocationFlag) ALLOCATE(IndexCauchy(BCFld%SizeCauchy))

        IF (process.EQ.0) THEN
            DO j=1,BCFld%SizeCauchy
                READ(u, '(I12,F15.10)')IndexCauchy(j),ValR
                CALL VecSetValue(BCFld%Cauchy(i),j-1,ValR,INSERT_VALUES,ierr)
            END DO
        END IF

        CALL VecAssemblyBegin(BCFld%Cauchy(i),ierr)
        CALL VecAssemblyEnd(BCFld%Cauchy(i),ierr)

    END DO

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetBC Stage] Boundary Conditions were satisfactorily created\n",ierr)

END SUBROUTINE GetBC_1

SUBROUTINE DestroyBC(BCFld,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(BoundaryConditions),INTENT(INOUT)  :: BCFld

    PetscInt                                :: i
    PetscLogStage                           :: Stage
    PetscBool                               :: Verbose

    CALL PetscLogStageRegister("DestroyBC", stage,ierr)
    CALL PetscLogStagePush(Stage,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[DestroyBC Stage] Inizialited\n",ierr)

    CALL VecDestroy(BCFld%Dirich,ierr)
    CALL VecDestroy(BCFld%Neumman,ierr)
    CALL VecDestroy(BCFld%Cauchy,ierr)
    DO i=1,BCFld%SizeTimeZone
        DEALLOCATE(BCFld%TimeZone(i)%Time)
    END DO
    DEALLOCATE(BCFld%TimeZone)

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[DestroyBC Stage] Finalized\n",ierr)

END SUBROUTINE DestroyBC

END MODULE ANISOFLOW_BoundaryConditions