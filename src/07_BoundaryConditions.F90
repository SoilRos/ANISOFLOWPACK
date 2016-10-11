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
    PetscLogDouble                          :: EventFlops=0.d0

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
    ELSEIF (InputType%BC.EQ.2) THEN
        CALL GetBC_2(Gmtry,BCFld,ierr)
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
    PetscInt                                :: u,i,j,ValI,CountTimeZone,CountDirich,CountSource,CountCauchy
    PetscInt                                :: LocalSizeDirich,LocalSizeSource,LocalSizeCauchy
    PetscInt,ALLOCATABLE                    :: IndexDirich(:),IndexSource(:),IndexCauchy(:)
    PetscReal                               :: ValR1,ValR2,DT
    CHARACTER(LEN=200)                      :: InputFileBC,InputDir,Route
    CHARACTER(LEN=200)                      :: aux
    PetscBool                               :: Verbose
    AO                                      :: AppOrd

    PARAMETER(u=01)

    CALL GetVerbose(Verbose,ierr)

    CALL GetInputDir(InputDir,ierr)
    CALL GetInputFileBC(InputFileBC,ierr)
    
    CALL DMDAGetAO(Gmtry%DataMngr,AppOrd,ierr)

    CALL MPI_Comm_rank(MPI_COMM_WORLD,process,ierr)

    IF (process.EQ.0) THEN
        Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFileBC))
        OPEN(u,FILE=TRIM(Route),STATUS='OLD',ACTION='READ')
        READ(u,*)aux,BCFld%SizeTimeZone
    END IF

    CALL MPI_Bcast(BCFld%SizeTimeZone,1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)
        
    ALLOCATE(BCFld%TimeZone(BCFld%SizeTimeZone))
    ALLOCATE(BCFld%SizeDirich(BCFld%SizeTimeZone))
    ALLOCATE(BCFld%Dirich(BCFld%SizeTimeZone))
    ALLOCATE(BCFld%SizeSource(BCFld%SizeTimeZone))
    ALLOCATE(BCFld%Source(BCFld%SizeTimeZone))
    ALLOCATE(BCFld%SizeCauchy(BCFld%SizeTimeZone))
    ALLOCATE(BCFld%CauchyC(BCFld%SizeTimeZone))
    ALLOCATE(BCFld%CauchyHe(BCFld%SizeTimeZone))
    ALLOCATE(BCFld%DirichIS(BCFld%SizeTimeZone))
    ALLOCATE(BCFld%SourceIS(BCFld%SizeTimeZone))
    ALLOCATE(BCFld%CauchyIS(BCFld%SizeTimeZone))

    DT=0.d0
    BCFld%SizeDirich(:)=0
    BCFld%SizeSource(:)=0
    BCFld%SizeCauchy(:)=0

    DO i=1,BCFld%SizeTimeZone
        ! TIME
        IF (process.EQ.0) THEN
            READ(u,*)aux,CountTimeZone,aux,BCFld%TimeZone(i)%SizeTime
        END IF

        CALL MPI_Bcast(BCFld%TimeZone(i)%SizeTime,1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)
        ALLOCATE(BCFld%TimeZone(i)%Time(BCFld%TimeZone(i)%SizeTime))

        IF (process.EQ.0) THEN
            READ(u,*)aux,DT
            DO j=1,BCFld%TimeZone(i)%SizeTime
                IF ((i.EQ.1).AND.(j.EQ.1)) THEN
                    BCFld%TimeZone(i)%Time(j)=DT
                ELSEIF (j.EQ.1) THEN
                    BCFld%TimeZone(i)%Time(j)=BCFld%TimeZone(i-1)%Time(BCFld%TimeZone(i-1)%SizeTime)+DT
                ELSE
                    BCFld%TimeZone(i)%Time(j)=BCFld%TimeZone(i)%Time(j-1)+DT
                END IF
            END DO
        END IF
        CALL MPI_Bcast(BCFld%TimeZone(i)%Time,BCFld%TimeZone(i)%SizeTime,MPI_DOUBLE, 0, PETSC_COMM_WORLD,ierr)

        ! Dirichlet
        IF (process.EQ.0) THEN
            READ(u,*)aux,CountDirich,aux,BCFld%SizeDirich(i)
        END IF

        LocalSizeDirich=BCFld%SizeDirich(i)
        CALL MPI_Bcast(BCFld%SizeDirich(i),1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)
        CALL MPI_Bcast(CountDirich,1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)

        IF (CountDirich.NE.i) THEN
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "[ERROR] Boundary condition file has to have a natural order on dirichlet identificators\n",ierr)
            STOP
            ! De hecho se pueden entrar pero es mejor restringirlo para que los archivos de entrada sean organizados
        END IF

        IF (BCFld%SizeDirich(i).GT.Gmtry%SizeTplgy(2)) THEN
!             CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
!             & "[ERROR] Boundary condition file doesn't has the same number of dirichlet entries as topology file\n",ierr)
!             STOP
            ! Si entra en este condicional significa que las condiciones de frontera no son permanentes.
        END IF

        CALL VecCreateMPI(PETSC_COMM_WORLD,LocalSizeDirich,BCFld%SizeDirich(i),BCFld%Dirich(CountDirich),ierr)
        
        ALLOCATE(IndexDirich(BCFld%SizeDirich(i)))

        IF (process.EQ.0) THEN
            DO j=1,BCFld%SizeDirich(i)
                READ(u, *)ValI,ValR1
                IndexDirich(j)=ValI-1
                ! Dirich is saved in its negative form
                CALL VecSetValue(BCFld%Dirich(CountDirich),j-1,-ValR1,INSERT_VALUES,ierr)
            END DO
        END IF
        
        CALL MPI_Bcast(IndexDirich,BCFld%SizeDirich(i),MPI_INT,0, PETSC_COMM_WORLD,ierr)
        CALL ISCreateGeneral(MPI_COMM_WORLD,BCFld%SizeDirich(i),IndexDirich,PETSC_COPY_VALUES,BCFld%DirichIS(CountDirich),ierr)

        DEALLOCATE(IndexDirich)

        CALL AOApplicationToPetscIS(AppOrd,BCFld%DirichIS(CountDirich),ierr)
        CALL VecAssemblyBegin(BCFld%Dirich(CountDirich),ierr)
        CALL VecAssemblyEnd(BCFld%Dirich(CountDirich),ierr)

        ! Source
        IF (process.EQ.0) THEN
            READ(u,*)aux,CountSource,aux,BCFld%SizeSource(i)
        END IF

        LocalSizeSource=BCFld%SizeSource(i)
        CALL MPI_Bcast(BCFld%SizeSource(i),1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)
        CALL MPI_Bcast(CountSource,1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)

        IF (CountSource.NE.i) THEN
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "[ERROR] Boundary condition file has to have a natural order on source identificators\n",ierr)
            STOP
            ! De hecho se pueden entrar pero es mejor restringirlo para que los archivos de entrada sean organizados
        END IF

        IF (BCFld%SizeSource(i).GT.Gmtry%SizeTplgy(3)) THEN
!             CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
!             & "[ERROR] Boundary condition file doesn't has the same number of dirichlet entries as topology file\n",ierr)
!             STOP
            ! Si entra en este condicional significa que las condiciones de frontera no son permanentes.
        END IF

        CALL VecCreateMPI(PETSC_COMM_WORLD,LocalSizeSource,BCFld%SizeSource(i),BCFld%Source(CountSource),ierr)

        ALLOCATE(IndexSource(BCFld%SizeSource(i)))

        IF (process.EQ.0) THEN
            DO j=1,BCFld%SizeSource(i)
                READ(u,*)ValI,ValR1
                IndexSource(j)=ValI-1
                CALL VecSetValue(BCFld%Source(CountSource),j-1,-ValR1,INSERT_VALUES,ierr)
            END DO
        END IF

        CALL MPI_Bcast(IndexSource,BCFld%SizeSource(i),MPI_INT,0, PETSC_COMM_WORLD,ierr)
        CALL ISCreateGeneral(MPI_COMM_WORLD,BCFld%SizeSource(i),IndexSource,PETSC_COPY_VALUES,BCFld%SourceIS(CountSource),ierr)

        DEALLOCATE(IndexSource)

        CALL AOApplicationToPetscIS(AppOrd,BCFld%SourceIS(CountSource),ierr)
        CALL VecAssemblyBegin(BCFld%Source(CountSource),ierr)
        CALL VecAssemblyEnd(BCFld%Source(CountSource),ierr)

        ! Cauchy
        IF (process.EQ.0) THEN
            READ(u,*)aux,CountCauchy,aux,BCFld%SizeCauchy(i)
        END IF

        LocalSizeCauchy=BCFld%SizeCauchy(i)
        CALL MPI_Bcast(BCFld%SizeCauchy(i),1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)
        CALL MPI_Bcast(CountCauchy,1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)

        IF (CountCauchy.NE.i) THEN
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "[ERROR] Boundary condition file has to have a natural order on cauchy identificators\n",ierr)
            STOP
            ! De hecho se pueden entrar pero es mejor restringirlo para que los archivos de entrada sean organizados
        END IF

        IF (BCFld%SizeCauchy(i).GT.Gmtry%SizeTplgy(4)) THEN
!             CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
!             & "[ERROR] Boundary condition file doesn't has the same number of dirichlet entries as topology file\n",ierr)
!             STOP
            ! Si entra en este condicional significa que las condiciones de frontera no son permanentes.
        END IF

        CALL VecCreateMPI(PETSC_COMM_WORLD,LocalSizeCauchy,BCFld%SizeCauchy(i),BCFld%CauchyC(CountCauchy),ierr)
        CALL VecDuplicate(BCFld%CauchyC(CountCauchy),BCFld%CauchyHe(CountCauchy),ierr)

        ALLOCATE(IndexCauchy(BCFld%SizeCauchy(i)))

        IF (process.EQ.0) THEN
            DO j=1,BCFld%SizeCauchy(i)
                READ(u,*)ValI,ValR1,ValR2
                IndexCauchy(j)=ValI-1
                ! Cauchy "He" is saved in its negative form
                CALL VecSetValue(BCFld%CauchyC(CountCauchy),j-1,ValR1,INSERT_VALUES,ierr)
                CALL VecSetValue(BCFld%CauchyHe(CountCauchy),j-1,-ValR2,INSERT_VALUES,ierr)
            END DO
        END IF

        CALL MPI_Bcast(IndexCauchy,BCFld%SizeCauchy(i),MPI_INT,0, PETSC_COMM_WORLD,ierr)
        CALL ISCreateGeneral(MPI_COMM_WORLD,BCFld%SizeCauchy(i),IndexCauchy,PETSC_COPY_VALUES,BCFld%CauchyIS(CountCauchy),ierr)

        DEALLOCATE(IndexCauchy)

        CALL AOApplicationToPetscIS(AppOrd,BCFld%CauchyIS(CountCauchy),ierr)
        CALL VecAssemblyBegin(BCFld%CauchyC(CountCauchy),ierr)
        CALL VecAssemblyEnd(BCFld%CauchyC(CountCauchy),ierr)
        CALL VecAssemblyBegin(BCFld%CauchyHe(CountCauchy),ierr)
        CALL VecAssemblyEnd(BCFld%CauchyHe(CountCauchy),ierr)

    END DO

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetBC Event] Boundary Conditions were satisfactorily created\n",ierr)

END SUBROUTINE GetBC_1

SUBROUTINE GetBC_2(Gmtry,BCFld,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,BoundaryConditions
    USE ANISOFLOW_Interface, ONLY : GetVerbose,GetInputDir,GetInputFileBC

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Geometry),INTENT(IN)               :: Gmtry
    TYPE(BoundaryConditions),INTENT(OUT)    :: BCFld

    PetscMPIInt                             :: process
    PetscInt                                :: u,i,j,ValI1,ValI2,ValI3,CountTimeZone,CountDirich,CountSource,CountCauchy,widthG(3)
    PetscInt                                :: LocalSizeDirich,LocalSizeSource,LocalSizeCauchy
    PetscInt,ALLOCATABLE                    :: IndexDirich(:),IndexSource(:),IndexCauchy(:)
    PetscReal                               :: ValR1,ValR2,DT
    CHARACTER(LEN=200)                      :: InputFileBC,InputDir,Route
    CHARACTER(LEN=200)                      :: aux
    PetscBool                               :: Verbose
    AO                                      :: AppOrd

    PARAMETER(u=01)

    CALL GetVerbose(Verbose,ierr)

    CALL GetInputDir(InputDir,ierr)
    CALL GetInputFileBC(InputFileBC,ierr)
    
    CALL DMDAGetAO(Gmtry%DataMngr,AppOrd,ierr)

    CALL MPI_Comm_rank(MPI_COMM_WORLD,process,ierr)

    IF (process.EQ.0) THEN
        Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFileBC))
        OPEN(u,FILE=TRIM(Route),STATUS='OLD',ACTION='READ')
        READ(u,*)aux,BCFld%SizeTimeZone
    END IF

    CALL MPI_Bcast(BCFld%SizeTimeZone,1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)
        
    ! It gets the global size from the geometry data manager.
    CALL DMDAGetInfo(Gmtry%DataMngr,PETSC_NULL_INTEGER,widthG(1),widthG(2),&
        & widthG(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,ierr)

    ALLOCATE(BCFld%TimeZone(BCFld%SizeTimeZone))
    ALLOCATE(BCFld%SizeDirich(BCFld%SizeTimeZone))
    ALLOCATE(BCFld%Dirich(BCFld%SizeTimeZone))
    ALLOCATE(BCFld%SizeSource(BCFld%SizeTimeZone))
    ALLOCATE(BCFld%Source(BCFld%SizeTimeZone))
    ALLOCATE(BCFld%SizeCauchy(BCFld%SizeTimeZone))
    ALLOCATE(BCFld%CauchyC(BCFld%SizeTimeZone))
    ALLOCATE(BCFld%CauchyHe(BCFld%SizeTimeZone))
    ALLOCATE(BCFld%DirichIS(BCFld%SizeTimeZone))
    ALLOCATE(BCFld%SourceIS(BCFld%SizeTimeZone))
    ALLOCATE(BCFld%CauchyIS(BCFld%SizeTimeZone))

    BCFld%SizeDirich(:)=0
    BCFld%SizeSource(:)=0
    BCFld%SizeCauchy(:)=0

    DO i=1,BCFld%SizeTimeZone
        ! TIME
        IF (process.EQ.0) THEN
            READ(u,*)aux,CountTimeZone,aux,BCFld%TimeZone(i)%SizeTime
        END IF

        CALL MPI_Bcast(BCFld%TimeZone(i)%SizeTime,1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)
        ALLOCATE(BCFld%TimeZone(i)%Time(BCFld%TimeZone(i)%SizeTime))

        IF (process.EQ.0) THEN
            READ(u,*)aux,DT
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
            READ(u,*)aux,CountDirich,aux,BCFld%SizeDirich(i)
        END IF
        
        LocalSizeDirich=BCFld%SizeDirich(i)
        CALL MPI_Bcast(BCFld%SizeDirich(i),1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)
        CALL MPI_Bcast(CountDirich,1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)

        IF (CountDirich.NE.i) THEN
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "[ERROR] Boundary condition file has to have a natural order on dirichlet identificators\n",ierr)
            STOP
            ! De hecho se pueden entrar pero es mejor restringirlo para que los archivos de entrada sean organizados
        END IF

        IF (BCFld%SizeDirich(i).GT.Gmtry%SizeTplgy(2)) THEN
!             CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
!             & "[ERROR] Boundary condition file doesn't has the same number of dirichlet entries as topology file\n",ierr)
!             STOP
            ! Si entra en este condicional significa que las condiciones de frontera no son permanentes.
        END IF

        CALL VecCreateMPI(PETSC_COMM_WORLD,LocalSizeDirich,BCFld%SizeDirich(i),BCFld%Dirich(CountDirich),ierr)
        
        ALLOCATE(IndexDirich(BCFld%SizeDirich(i)))

        IF (process.EQ.0) THEN
            DO j=1,BCFld%SizeDirich(i)
                READ(u,*)ValI1,ValI2,ValI3,ValR1
                ValI1=widthG(1)*widthG(2)*(ValI3-1)+widthG(1)*(ValI2-1)+ValI1
                IndexDirich(j)=ValI1-1
                ! Dirich is saved in its negative form
                CALL VecSetValue(BCFld%Dirich(CountDirich),j-1,-ValR1,INSERT_VALUES,ierr)
            END DO
        END IF
        
        CALL MPI_Bcast(IndexDirich,BCFld%SizeDirich(i),MPI_INT,0, PETSC_COMM_WORLD,ierr)
        CALL ISCreateGeneral(MPI_COMM_WORLD,BCFld%SizeDirich(i),IndexDirich,PETSC_COPY_VALUES,BCFld%DirichIS(CountDirich),ierr)

        DEALLOCATE(IndexDirich)

        CALL AOApplicationToPetscIS(AppOrd,BCFld%DirichIS(CountDirich),ierr)
        CALL VecAssemblyBegin(BCFld%Dirich(CountDirich),ierr)
        CALL VecAssemblyEnd(BCFld%Dirich(CountDirich),ierr)
        
        ! Source
        IF (process.EQ.0) THEN
            READ(u,*)aux,CountSource,aux,BCFld%SizeSource(i)
        END IF

        LocalSizeSource=BCFld%SizeSource(i)
        CALL MPI_Bcast(BCFld%SizeSource(i),1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)
        CALL MPI_Bcast(CountSource,1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)

        IF (CountSource.NE.i) THEN
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "[ERROR] Boundary condition file has to have a natural order on source identificators\n",ierr)
            STOP
            ! De hecho se pueden entrar pero es mejor restringirlo para que los archivos de entrada sean organizados
        END IF

        IF (BCFld%SizeSource(i).GT.Gmtry%SizeTplgy(3)) THEN
!             CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
!             & "[ERROR] Boundary condition file doesn't has the same number of dirichlet entries as topology file\n",ierr)
!             STOP
            ! Si entra en este condicional significa que las condiciones de frontera no son permanentes.
        END IF

        CALL VecCreateMPI(PETSC_COMM_WORLD,LocalSizeSource,BCFld%SizeSource(i),BCFld%Source(CountSource),ierr)

        ALLOCATE(IndexSource(BCFld%SizeSource(i)))

        IF (process.EQ.0) THEN
            DO j=1,BCFld%SizeSource(i)
                READ(u,*)ValI1,ValI2,ValI3,ValR1
                ValI1=widthG(1)*widthG(2)*(ValI3-1)+widthG(1)*(ValI2-1)+ValI1
                IndexSource(j)=ValI1-1
                CALL VecSetValue(BCFld%Source(CountSource),j-1,-ValR1,INSERT_VALUES,ierr)
            END DO
        END IF

        CALL MPI_Bcast(IndexSource,BCFld%SizeSource(i),MPI_INT,0, PETSC_COMM_WORLD,ierr)
        CALL ISCreateGeneral(MPI_COMM_WORLD,BCFld%SizeSource(i),IndexSource,PETSC_COPY_VALUES,BCFld%SourceIS(CountSource),ierr)

        DEALLOCATE(IndexSource)

        CALL AOApplicationToPetscIS(AppOrd,BCFld%SourceIS(CountSource),ierr)
        CALL VecAssemblyBegin(BCFld%Source(CountSource),ierr)
        CALL VecAssemblyEnd(BCFld%Source(CountSource),ierr)

        ! Cauchy
        IF (process.EQ.0) THEN
            READ(u,*)aux,CountCauchy,aux,BCFld%SizeCauchy(i)
        END IF

        LocalSizeCauchy=BCFld%SizeCauchy(i)
        CALL MPI_Bcast(BCFld%SizeCauchy(i),1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)
        CALL MPI_Bcast(CountCauchy,1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)

        IF (CountCauchy.NE.i) THEN
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "[ERROR] Boundary condition file has to have a natural order on cauchy identificators\n",ierr)
            STOP
            ! De hecho se pueden entrar pero es mejor restringirlo para que los archivos de entrada sean organizados
        END IF

        IF (BCFld%SizeCauchy(i).GT.Gmtry%SizeTplgy(4)) THEN
!             CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
!             & "[ERROR] Boundary condition file doesn't has the same number of dirichlet entries as topology file\n",ierr)
!             STOP
            ! Si entra en este condicional significa que las condiciones de frontera no son permanentes.
        END IF

        CALL VecCreateMPI(PETSC_COMM_WORLD,LocalSizeCauchy,BCFld%SizeCauchy(i),BCFld%CauchyC(CountCauchy),ierr)
        CALL VecDuplicate(BCFld%CauchyC(CountCauchy),BCFld%CauchyHe(CountCauchy),ierr)

        ALLOCATE(IndexCauchy(BCFld%SizeCauchy(i)))

        IF (process.EQ.0) THEN
            DO j=1,BCFld%SizeCauchy(i)
                READ(u,*)ValI1,ValI2,ValI3,ValR1,ValR2
                ValI1=widthG(1)*widthG(2)*(ValI3-1)+widthG(1)*(ValI2-1)+ValI1
                IndexCauchy(j)=ValI1-1
                ! Cauchy "He" is saved in its negative form
                CALL VecSetValue(BCFld%CauchyC(CountCauchy),j-1,ValR1,INSERT_VALUES,ierr)
                CALL VecSetValue(BCFld%CauchyHe(CountCauchy),j-1,-ValR2,INSERT_VALUES,ierr)
            END DO
        END IF

        CALL MPI_Bcast(IndexCauchy,BCFld%SizeCauchy(i),MPI_INT,0, PETSC_COMM_WORLD,ierr)
        CALL ISCreateGeneral(MPI_COMM_WORLD,BCFld%SizeCauchy(i),IndexCauchy,PETSC_COPY_VALUES,BCFld%CauchyIS(CountCauchy),ierr)

        DEALLOCATE(IndexCauchy)

        CALL AOApplicationToPetscIS(AppOrd,BCFld%CauchyIS(CountCauchy),ierr)
        CALL VecAssemblyBegin(BCFld%CauchyC(CountCauchy),ierr)
        CALL VecAssemblyEnd(BCFld%CauchyC(CountCauchy),ierr)
        CALL VecAssemblyBegin(BCFld%CauchyHe(CountCauchy),ierr)
        CALL VecAssemblyEnd(BCFld%CauchyHe(CountCauchy),ierr)

    END DO

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetBC Event] Boundary Conditions were satisfactorily created\n",ierr)

END SUBROUTINE GetBC_2

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
    PetscLogDouble                          :: EventFlops=0.d0

    ClassName="BC"
    CALL PetscClassIdRegister(ClassName,ClassID,ierr)
    EventName="DestroyBC"
    CALL PetscLogEventRegister(EventName,ClassID,Event,ierr)
    CALL PetscLogEventBegin(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Inizialited\n",ierr)

    CALL VecDestroy(BCFld%Dirich,ierr)
    CALL VecDestroy(BCFld%Source,ierr)
    CALL VecDestroy(BCFld%CauchyC,ierr)
    CALL VecDestroy(BCFld%CauchyHe,ierr)
    DO i=1,BCFld%SizeTimeZone
        DEALLOCATE(BCFld%TimeZone(i)%Time)
    END DO
    DEALLOCATE(BCFld%TimeZone)

    CALL ISDestroy(BCFld%DirichIS,ierr)
    CALL ISDestroy(BCFld%SourceIS,ierr)
    CALL ISDestroy(BCFld%CauchyIS,ierr)

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Finalized\n",ierr)
    
    CALL PetscLogFlops(EventFlops,ierr)
    CALL PetscLogEventEnd(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

END SUBROUTINE DestroyBC

END MODULE ANISOFLOW_BoundaryConditions