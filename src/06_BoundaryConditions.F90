MODULE ANISOFLOW_BoundaryConditions

    ! USE ANISOFLOW_Interface, ONLY : GetVerbose
    ! USE ANISOFLOW_Types, ONLY : Geometry,BoundaryConditions

    IMPLICIT NONE

CONTAINS

SUBROUTINE GetBC(Gmtry,BCFld,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,BoundaryConditions,InputTypeVar
    USE ANISOFLOW_Interface, ONLY : GetVerbose,GetInputType

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Geometry),INTENT(IN)               :: Gmtry
    TYPE(BoundaryConditions),INTENT(OUT)    :: BCFld

    TYPE(InputTypeVar)                      :: InputType
    CHARACTER(LEN=200)                      :: EventName,ClassName
    PetscBool                               :: Verbose
    PetscLogEvent                           :: Event
    PetscClassId                            :: ClassID
    PetscLogDouble                          :: EventFlops

    ClassName="BC"
    CALL PetscClassIdRegister(ClassName,ClassID,ierr)
    EventName="GetBC"
    CALL PetscLogEventRegister(EventName,ClassID,Event,ierr)
    CALL PetscLogEventBegin(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Inizialited\n",ierr)
    
    CALL GetInputType(InputType,ierr)

    IF (InputType%BC.EQ.1) THEN
        CALL GetBC_1(Gmtry,BCFld,ierr)
    ELSE 
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "[ERROR] Boundary Condition InputType wrong\n",ierr)
        STOP
    END IF

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Finalized\n",ierr)
    
    CALL PetscLogFlops(EventFlops,ierr)
    CALL PetscLogEventEnd(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)


END SUBROUTINE GetBC

SUBROUTINE GetBC_1(Gmtry,BCFld,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,BoundaryConditions
    USE ANISOFLOW_Interface, ONLY : GetVerbose,GetInputDir,GetInputFileBC

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Geometry),INTENT(IN)               :: Gmtry
    TYPE(BoundaryConditions),INTENT(OUT)    :: BCFld

    PetscMPIInt                             :: process
    PetscInt                                :: u,i,j,CountTimeZone,CountDirich,CountSource,CountCauchy
    PetscInt                                :: SizeDirichIS,SizeSourceIS,SizeCauchyIS
    PetscInt,ALLOCATABLE                    :: IndexDirich(:),IndexSource(:),IndexCauchy(:)
    PetscReal                               :: ValR,DT
    CHARACTER(LEN=200)                      :: InputFileBC,InputDir,Route
    CHARACTER(LEN=20)                       :: TextTimeZones
    CHARACTER(LEN=11)                       :: Text1TimeZone
    CHARACTER(LEN=3)                        :: Text2TimeZone,Text2Dirich,Text2Source,Text2Cauchy
    CHARACTER(LEN=4)                        :: TextDT
    CHARACTER(LEN=27)                       :: Text1Dirich
    CHARACTER(LEN=25)                       :: Text1Source
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
    ALLOCATE(BCFld%Source(BCFld%SizeTimeZone))
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
            DO j=1,BCFld%TimeZone(i)%SizeTime
                IF ((i.EQ.1).AND.(j.EQ.1)) THEN
                    BCFld%TimeZone(1)%Time(1)=0.0
                ELSEIF (j.EQ.1) THEN
                    BCFld%TimeZone(i)%Time(1)=BCFld%TimeZone(i-1)%Time(BCFld%TimeZone(i-1)%SizeTime)+DT
                ELSE
                    BCFld%TimeZone(i)%Time(j)=BCFld%TimeZone(i)%Time(j-1)+DT
                END IF
            END DO
        END IF
        CALL MPI_Bcast(BCFld%TimeZone(i)%Time,BCFld%TimeZone(i)%SizeTime,MPI_DOUBLE, 0, PETSC_COMM_WORLD,ierr)

        ! Dirichlet
        IF (process.EQ.0) THEN
            READ(u,'(A27,I8,A3,I8)')Text1Dirich,CountDirich,Text2Dirich,BCFld%SizeDirich
            ! Imprimir errores de lectura de texto
        END IF
        
        CALL MPI_Bcast(BCFld%SizeDirich,1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)
        CALL ISGetLocalSize(Gmtry%DirichIS,SizeDirichIS,ierr)
        
        IF (BCFld%SizeDirich.NE.SizeDirichIS) THEN
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "[ERROR] Boundary condition file doesn't has the same number of dirichlet entries as topology file\n",ierr)
            STOP
        END IF

        CALL VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,BCFld%SizeDirich,BCFld%Dirich(i),ierr)
        
        AllocationFlag=ALLOCATED(IndexDirich)
        IF (.NOT.AllocationFlag) ALLOCATE(IndexDirich(BCFld%SizeDirich))

        IF (process.EQ.0) THEN
            DO j=1,BCFld%SizeDirich
                READ(u, '(I12,F15.10)')IndexDirich(j),ValR
                CALL VecSetValue(BCFld%Dirich(i),j-1,-ValR,INSERT_VALUES,ierr)
            END DO
        END IF

        CALL VecAssemblyBegin(BCFld%Dirich(i),ierr)
        CALL VecAssemblyEnd(BCFld%Dirich(i),ierr)
        
        ! Source
        IF (process.EQ.0) THEN
            READ(u,'(A25,I8,A3,I8)')Text1Source,CountSource,Text2Source,BCFld%SizeSource
            ! Imprimir errores de lectura de texto
        END IF

        CALL MPI_Bcast(BCFld%SizeSource,1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)
        CALL ISGetLocalSize(Gmtry%SourceIS,SizeSourceIS,ierr)

        IF (BCFld%SizeSource.NE.SizeSourceIS) THEN
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "[ERROR] Boundary condition file doesn't has the same number of source entries as topology file\n",ierr)
            STOP
        END IF

        CALL VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,BCFld%SizeSource,BCFld%Source(i),ierr)

        AllocationFlag=ALLOCATED(IndexSource)
        IF (.NOT.AllocationFlag) ALLOCATE(IndexSource(BCFld%SizeSource))

        IF (process.EQ.0) THEN
            DO j=1,BCFld%SizeSource
                READ(u, '(I12,F15.10)')IndexSource(j),ValR
                CALL VecSetValue(BCFld%Source(i),j-1,-ValR,INSERT_VALUES,ierr)
            END DO
        END IF

        CALL VecAssemblyBegin(BCFld%Source(i),ierr)
        CALL VecAssemblyEnd(BCFld%Source(i),ierr)

        ! Cauchy
        IF (process.EQ.0) THEN
            READ(u,'(A24,I8,A3,I8)')Text1Cauchy,CountCauchy,Text2Cauchy,BCFld%SizeCauchy
            ! Imprimir errores de lectura de texto
        END IF

        CALL MPI_Bcast(BCFld%SizeCauchy,1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)
        CALL ISGetLocalSize(Gmtry%CauchyIS,SizeCauchyIS,ierr)

        IF (BCFld%SizeCauchy.NE.SizeCauchyIS) THEN
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "[ERROR] Boundary condition file doesn't has the same number of cauchy entries as topology file\n",ierr)
            STOP
        END IF
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

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetBC Event] Boundary Conditions were satisfactorily created\n",ierr)

END SUBROUTINE GetBC_1

SUBROUTINE DestroyBC(BCFld,ierr)

    USE ANISOFLOW_Types, ONLY : BoundaryConditions
    USE ANISOFLOW_Interface, ONLY : GetVerbose

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(BoundaryConditions),INTENT(INOUT)  :: BCFld

    PetscInt                                :: i
    PetscBool                               :: Verbose
    CHARACTER(LEN=200)                      :: EventName,ClassName
    PetscLogEvent                           :: Event
    PetscClassId                            :: ClassID
    PetscLogDouble                          :: EventFlops

    ClassName="BC"
    CALL PetscClassIdRegister(ClassName,ClassID,ierr)
    EventName="DestroyBC"
    CALL PetscLogEventRegister(EventName,ClassID,Event,ierr)
    CALL PetscLogEventBegin(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Inizialited\n",ierr)

    CALL VecDestroy(BCFld%Dirich,ierr)
    CALL VecDestroy(BCFld%Source,ierr)
    CALL VecDestroy(BCFld%Cauchy,ierr)
    DO i=1,BCFld%SizeTimeZone
        DEALLOCATE(BCFld%TimeZone(i)%Time)
    END DO
    DEALLOCATE(BCFld%TimeZone)

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Finalized\n",ierr)
    
    CALL PetscLogFlops(EventFlops,ierr)
    CALL PetscLogEventEnd(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

END SUBROUTINE DestroyBC

END MODULE ANISOFLOW_BoundaryConditions