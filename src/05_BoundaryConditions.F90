MODULE ANISOFLOW_BoundaryConditions

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
    TYPE(BoundaryConditions),INTENT(INOUT)  :: BCFld

    ! CHARACTER(LEN=200)                      :: InputFileBC,InputDir,Route
    TYPE(InputTypeVar)                      :: InputType
    ! PetscMPIInt                             :: process
    ! PetscInt                                :: u,i,DirichLen
    ! PetscReal                               :: ValR
    ! PetscInt,ALLOCATABLE                    :: IndexDirich(:)
    ! IS                                      :: TempIS

    CALL GetInputType(InputType,ierr)

    IF (InputType%BC.EQ.1) THEN
        CALL GetBC_1(Gmtry,BCFld,ierr)
    ELSE 
        PRINT*,"ERROR"
    END IF

END SUBROUTINE GetBC

SUBROUTINE GetBC_1(Gmtry,BCFld,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Geometry),INTENT(IN)               :: Gmtry
    TYPE(BoundaryConditions),INTENT(INOUT)  :: BCFld

    PetscMPIInt                             :: process
    PetscInt                                :: u,i,SizeTimeZones,SizeTimeZone
    CHARACTER(LEN=200)                      :: InputFileBC,InputDir,Route
    CHARACTER(LEN=11)                       :: TextTimeZones,TextTimeZone

    PARAMETER(u=01)

    CALL GetInputDir(InputDir,ierr)
    CALL GetInputType(InputType,ierr)
    CALL GetInputFileBC(InputFileBC,ierr)
    
    CALL MPI_Comm_rank(MPI_COMM_WORLD,process,ierr)

    IF (process.EQ.0) THEN
        Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFileBC))
        OPEN(u,FILE=TRIM(Route),STATUS='OLD',ACTION='READ')
        READ(u,'(A11,I10)')TextTimeZones,SizeTimeZones
        IF (ADJUSTL(TRIM(TextTimeZones).NE."Timezones:") THEN
            !PRINT ERROR
        END IF
        ALLOCATE(BCFld%TimeZone(SizeTimeZones))
        DO i=1,SizeTimeZones
            READ(u,'(A11,I10)')TextTimeZone,SizeTimeZone
        END DO
        READ(u, '(I10)')

            ALLOCATE(IndexDirich(DirichLen))
        END IF


END SUBROUTINE GetBC_1









        CALL MPI_Comm_rank(MPI_COMM_WORLD,process,ierr)
        ! DirichLen=0
        IF (process.EQ.0) THEN
            Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFileBC))
            OPEN(u,FILE=TRIM(Route),STATUS='OLD',ACTION='READ')
            READ(u,*)"Timezones: "
            READ(u, '(I10)')DirichLen

            ALLOCATE(IndexDirich(DirichLen))
        END IF

        CALL MPI_Bcast(DirichLen,1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)
        CALL VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,DirichLen,             &
            & BCFld%Dirich(1),ierr)
        IF (process.EQ.0) THEN
            DO i=1,DirichLen
                READ(u, '(I12,F15.10)')IndexDirich(i),ValR
                CALL VecSetValue(BCFld%Dirich(1),i-1,ValR,INSERT_VALUES,  &
                    & ierr)
            END DO
            CLOSE(u)
        ELSE
            ALLOCATE(IndexDirich(DirichLen))
        END IF

        CALL VecAssemblyBegin(BCFld%Dirich(1),ierr)
        CALL VecAssemblyEnd(BCFld%Dirich(1),ierr)

        CALL MPI_Bcast(DirichLen,1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)
        CALL MPI_Bcast(IndexDirich,DirichLen,MPI_INT, 0, PETSC_COMM_WORLD,ierr)

        CALL ISCreateGeneral(PETSC_COMM_WORLD,DirichLen,IndexDirich,            &
            & PETSC_COPY_VALUES,TempIS,ierr)

        ! CALL CheckDirichIS() ! It must compare topologic and BC IS
    END IF




END SUBROUTINE GetSteadyBCDirichlet

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Geometry),INTENT(IN)               :: Gmtry
    TYPE(BoundaryConditions),INTENT(INOUT)  :: BCFld




SUBROUTINE GetBC(Gmtry,BCFld,ierr)
    
    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Geometry),INTENT(IN)               :: Gmtry
    TYPE(BoundaryConditions),INTENT(OUT)    :: BCFld

    TYPE(RunOptionsVar)                     :: RunOptions

    PetscReal                               :: Value
    PetscInt                                :: i

    CALL GetRunOptions(RunOptions,ierr)

    ALLOCATE(BCFld%TimeZone(1))                                         ! TEST!!!!
    ALLOCATE(BCFld%TimeZone(1)%Time(10))
    DO i=0,9                                                        ! TEST!!!!
        BCFld%TimeZone(1)%Time(i+1)=REAL(i)
    END DO                                                          ! TEST!!!!

    ! IF (RunOptions%Time) THEN
        ! CALL GetTransientBC(Gmtry,BCFld,ierr)
    ! ELSE
        CALL GetSteadyBC(Gmtry,BCFld,ierr)
    ! END IF

END SUBROUTINE GetBC

SUBROUTINE GetSteadyBC(Gmtry,BCFld,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Geometry),INTENT(IN)               :: Gmtry
    TYPE(BoundaryConditions),INTENT(OUT)    :: BCFld

    ALLOCATE(BCFld%Dirich(1))
    CALL GetSteadyBCDirichlet(Gmtry,BCFld,ierr)
    ! CALL GetSteadyBCCauchy(Gmtry,BCFld,ierr)

END SUBROUTINE GetSteadyBC

SUBROUTINE GetSteadyBCDirichlet(Gmtry,BCFld,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Geometry),INTENT(IN)               :: Gmtry
    TYPE(BoundaryConditions),INTENT(INOUT)  :: BCFld

    CHARACTER(LEN=200)                      :: InputFileBC,InputDir,Route
    TYPE(InputTypeVar)                      :: InputType
    PetscMPIInt                             :: process
    PetscInt                                :: u,i,DirichLen
    PetscReal                               :: ValR
    PetscInt,ALLOCATABLE                    :: IndexDirich(:)
    IS                                      :: TempIS

    PARAMETER(u=01)
    ValR=0
    CALL GetInputDir(InputDir,ierr)
    CALL GetInputType(InputType,ierr)
    CALL GetInputFileBC(InputFileBC,ierr)

    IF (InputType%BC.EQ.1) THEN

        CALL MPI_Comm_rank(MPI_COMM_WORLD,process,ierr)
        DirichLen=0
        IF (process.EQ.0) THEN
            Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFileBC))
            OPEN(u,FILE=TRIM(Route),STATUS='OLD',ACTION='READ')
            READ(u,*)
            READ(u, '(I10)')DirichLen

            ALLOCATE(IndexDirich(DirichLen))
        END IF

        CALL MPI_Bcast(DirichLen,1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)
        CALL VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,DirichLen,             &
            & BCFld%Dirich(1),ierr)
        IF (process.EQ.0) THEN
            DO i=1,DirichLen
                READ(u, '(I12,F15.10)')IndexDirich(i),ValR
                CALL VecSetValue(BCFld%Dirich(1),i-1,ValR,INSERT_VALUES,  &
                    & ierr)
            END DO
            CLOSE(u)
        ELSE
            ALLOCATE(IndexDirich(DirichLen))
        END IF

        CALL VecAssemblyBegin(BCFld%Dirich(1),ierr)
        CALL VecAssemblyEnd(BCFld%Dirich(1),ierr)

        CALL MPI_Bcast(DirichLen,1,MPI_INT, 0, PETSC_COMM_WORLD,ierr)
        CALL MPI_Bcast(IndexDirich,DirichLen,MPI_INT, 0, PETSC_COMM_WORLD,ierr)

        CALL ISCreateGeneral(PETSC_COMM_WORLD,DirichLen,IndexDirich,            &
            & PETSC_COPY_VALUES,TempIS,ierr)

        ! CALL CheckDirichIS() ! It must compare topologic and BC IS
    END IF


END SUBROUTINE GetSteadyBCDirichlet

SUBROUTINE GetSteadyBCCauchy(Gmtry,BCFld,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Geometry),INTENT(IN)               :: Gmtry
    TYPE(BoundaryConditions),INTENT(INOUT)  :: BCFld

    CHARACTER(LEN=200)                      :: InputFileBC,InputDir,Route
    TYPE(InputTypeVar)                      :: InputType
    PetscInt                                :: u

    PARAMETER(u=01)

    CALL GetInputDir(InputDir,ierr)
    CALL GetInputType(InputType,ierr)
    CALL GetInputFileBC(InputFileBC,ierr)

    IF (InputType%BC.EQ.1) THEN

        ! InputType SeteadyBC equal 1 don't use cauchy BC

    END IF


END SUBROUTINE GetSteadyBCCauchy

SUBROUTINE BCDestroy(BCFld,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    TYPE(BoundaryConditions),INTENT(IN) :: BCFld ! INOUT or IN?

    CALL VecDestroy(BCFld%Dirich,ierr)
    ! CALL VecDestroy(BCFld%Cauchy,ierr)

END SUBROUTINE BCDestroy

END MODULE ANISOFLOW_BoundaryConditions