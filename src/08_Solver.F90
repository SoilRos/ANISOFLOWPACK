MODULE ANISOFLOW_Solver

    IMPLICIT NONE

CONTAINS

SUBROUTINE SolveSystem(Gmtry,PptFld,BCFld,A,b,x,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,PropertyField,BoundaryConditions,RunOptionsVar
    USE ANISOFLOW_Interface, ONLY : GetRunOptions
    USE ANISOFLOW_View, ONLY : ViewSolution
    USE ANISOFLOW_BuildSystem, ONLY : GetSystem

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscviewer.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Geometry),INTENT(IN)               :: Gmtry
    TYPE(PropertyField),INTENT(IN)          :: PptFld
    TYPE(BoundaryConditions),INTENT(IN)     :: BCFld
    Mat,INTENT(INOUT)                       :: A
    Vec,INTENT(INOUT)                       :: b,x

    KSP                                     :: Solver
    TYPE(RunOptionsVar)                     :: RunOptions
    Vec                                     :: diagA
    CHARACTER(LEN=200)                      :: Name,CharCount

    PetscInt                                :: i,j,initTime,Count
    PetscReal                               :: zero=0.0

    CALL KSPCreate(PETSC_COMM_WORLD,Solver,ierr)
    CALL GetRunOptions(RunOptions,ierr)

    IF (RunOptions%Time) THEN

        CALL VecDuplicate(x,diagA,ierr)
        CALL MatGetDiagonal(A,diagA,ierr)
        Count=102
        DO i=1,BCFld%SizeTimeZone

            initTime=1
            IF (i.EQ.1) initTime=2 



            DO j=initTime,BCFld%TimeZone(i)%SizeTime

                CALL GetSystem(Gmtry,PptFld,BCFld,i,j,A,b,x,ierr)
                CALL KSPSetOperators(Solver,A,A,ierr)
                CALL KSPSetTolerances(Solver,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,    &
                    & PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
                CALL KSPSetFromOptions(Solver,ierr)
                CALL KSPSolve(Solver,b,x,ierr)
            
                WRITE(CharCount,*)Count
                Name="ANISOFLOW_sol_"//TRIM(ADJUSTL(CharCount))//".h5"
                Name=ADJUSTL(Name)
                print*,Name
                CALL ViewSolution(x,Name,ierr)
                Count=Count+1

                CALL MatDiagonalSet(A,diagA,INSERT_VALUES,ierr)
                CALL VecSet(b,zero,ierr)
            END DO



        END DO

    ELSE

        CALL KSPSetOperators(Solver,A,A,ierr)
        CALL KSPSetTolerances(Solver,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,    &
            & PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
        CALL KSPSetFromOptions(Solver,ierr)
        CALL KSPSolve(Solver,b,x,ierr)
        Name="ANISOFLOW_sol.h5"
        CALL ViewSolution(x,Name,ierr)

    END IF



    CALL KSPDestroy(Solver,ierr)

END SUBROUTINE SolveSystem

END MODULE ANISOFLOW_Solver
