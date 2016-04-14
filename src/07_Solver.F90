MODULE ANISOFLOW_Solver

    IMPLICIT NONE

CONTAINS

SUBROUTINE SolveSystem(Gmtry,BCFld,A,b,x,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,BoundaryConditions
    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscksp.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Geometry),INTENT(IN)               :: Gmtry
    TYPE(BoundaryConditions),INTENT(IN)     :: BCFld
    Mat,INTENT(INOUT)                       :: A
    Vec,INTENT(INOUT)                       :: b,x

    KSP                                     :: Solver
    TYPE(RunOptionsVar)                     :: RunOptions

    CALL KSPCreate(PETSC_COMM_WORLD,Solver,ierr)
    CALL GetRunOptions(RunOptions,ierr)

    IF (RunOptions%Time) THEN

        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,             &
            & "ERROR: Transitory not implemented yet, please use steady simulation for now",ierr)
        STOP

    ELSE

        CALL KSPSetOperators(Solver,A,A,ierr)
        CALL KSPSetTolerances(Solver,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,    &
            & PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
        CALL KSPSetFromOptions(Solver,ierr)
        CALL KSPSolve(Solver,b,x,ierr)

    END IF



    CALL KSPDestroy(Solver,ierr)

END SUBROUTINE SolveSystem

END MODULE ANISOFLOW_Solver
