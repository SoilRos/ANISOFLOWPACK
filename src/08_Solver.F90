MODULE ANISOFLOW_Solver

    IMPLICIT NONE

CONTAINS

SUBROUTINE SolveSystem(Gmtry,PptFld,BCFld,A,b,x,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,PropertiesField,BoundaryConditions,RunOptionsVar
    USE ANISOFLOW_Interface, ONLY : GetRunOptions, GetVerbose
    USE ANISOFLOW_View, ONLY : ViewSolution
    USE ANISOFLOW_BuildSystem, ONLY : GetSystem

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscviewer.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Geometry),INTENT(INOUT)            :: Gmtry
    TYPE(PropertiesField),INTENT(IN)        :: PptFld
    TYPE(BoundaryConditions),INTENT(IN)     :: BCFld
    Mat,INTENT(INOUT)                       :: A
    Vec,INTENT(INOUT)                       :: b,x

    KSP                                     :: Solver
    Vec                                     :: diagA
    CHARACTER(LEN=200)                      :: Name,CharCount,Chari,Charj,CharMsg

    PetscInt                                :: i,j,Count
    PetscReal                               :: zero=0.0

    TYPE(RunOptionsVar)                     :: RunOptions
    CHARACTER(LEN=200)                      :: EventName,ClassName
    PetscBool                               :: Verbose
    PetscLogEvent                           :: Event
    PetscClassId                            :: ClassID
    PetscLogDouble                          :: EventFlops=0.d0


    ClassName="System"
    CALL PetscClassIdRegister(ClassName,ClassID,ierr)
    EventName="SolveSystem"
    CALL PetscLogEventRegister(EventName,ClassID,Event,ierr)
    CALL PetscLogEventBegin(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Inizialited\n",ierr)
    

    CALL KSPCreate(PETSC_COMM_WORLD,Solver,ierr)
    CALL GetRunOptions(RunOptions,ierr)

    IF (RunOptions%Time) THEN ! Transitory

        CALL VecDuplicate(x,diagA,ierr)
        CALL MatGetDiagonal(A,diagA,ierr)
        Count=1
        DO i=1,BCFld%SizeTimeZone
            DO j=1,BCFld%TimeZone(i)%SizeTime
                
                WRITE(Chari,*)i
                WRITE(Charj,*)j

                CharMsg="["//TRIM(ADJUSTL(Chari))//":"//TRIM(ADJUSTL(Charj))//"]"
!                 print*,CharMsg
                IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Transitory iteration "//TRIM(CharMsg)// " inizialited\n",ierr)
                
                IF (i.NE.1) THEN
                    CALL GetSystem(Gmtry,PptFld,BCFld,i,j,A,b,x,ierr)
                    CALL KSPSetOperators(Solver,A,A,ierr)
                    CALL KSPSetTolerances(Solver,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,    &
                        & PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
                    CALL KSPSetFromOptions(Solver,ierr)
                    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] PETSc solver monitor:\n",ierr)
                    CALL KSPSolve(Solver,b,x,ierr)
                END IF
                WRITE(CharCount,*)Count
                Name="ANISOFLOW_Sol_"//TRIM(ADJUSTL(CharCount))
                Name=ADJUSTL(Name)
                CALL ViewSolution(x,Name,EventName,ierr)
                Count=Count+1
                CALL MatDiagonalSet(A,diagA,INSERT_VALUES,ierr)
                CALL VecSet(b,zero,ierr)

                IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Transitory iteration "//TRIM(CharMsg)// " finalized\n",ierr)

                ! Re-create the solver just to close any viewer opened before.
                CALL KSPDestroy(Solver,ierr)
                CALL KSPCreate(PETSC_COMM_WORLD,Solver,ierr)

            END DO
        END DO
    ELSE ! Steady

        IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Steady solution inizialited\n",ierr)
   
        Name="ANISOFLOW_b"
        CALL ViewSolution(b,Name,EventName,ierr)

        CALL KSPSetOperators(Solver,A,A,ierr)
        CALL KSPSetTolerances(Solver,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,    &
            & PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
        CALL KSPSetFromOptions(Solver,ierr)
        CALL KSPSolve(Solver,b,x,ierr)

        Name="ANISOFLOW_Sol"
        CALL ViewSolution(x,Name,EventName,ierr)

        IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Steady solution finalized\n",ierr)

    END IF

    CALL KSPDestroy(Solver,ierr)

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Finalized\n",ierr)
    
    CALL PetscLogFlops(EventFlops,ierr)
    CALL PetscLogEventEnd(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)


END SUBROUTINE SolveSystem

END MODULE ANISOFLOW_Solver
