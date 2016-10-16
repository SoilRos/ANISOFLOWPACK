MODULE ANISOFLOW_Solver

    IMPLICIT NONE

CONTAINS

SUBROUTINE SolveSystem(Gmtry,PptFld,BCFld,A,b,x,ierr)

    USE ANISOFLOW_Types,       ONLY : Geometry,PropertiesField,      &
                                    & BoundaryConditions,RunOptionsVar
    USE ANISOFLOW_Interface,   ONLY : GetRunOptions, GetVerbose
    USE ANISOFLOW_Geometry,    ONLY : UpdateTplgy
    USE ANISOFLOW_View,        ONLY : ViewSolution
    USE ANISOFLOW_BuildSystem, ONLY : ApplyDirichlet,ApplySource,    &
                                    & ApplyCauchy,ApplyTimeDiff,     &
                                    & GetInitSol

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
    KSPType                                 :: SolverType
    KSPConvergedReason                      :: SolverConvergedReason
    PC                                      :: SolverPC
    PCType                                  :: SolverPCType
!     Vec                                     :: diagA
    Mat                                     :: copyA
    CHARACTER(LEN=200)                      :: ViewName,CharCount,   &
                                             & Chari,Charj,CharMsg,  &
                                             & CharTolR,CharTolAbs,  &
                                             & CharTolD,CharMaxIts,  &
                                             & CharSolverConvergedReason
    
    PetscInt                                :: i,j,Count,MaxIts
    PetscReal                               :: TolR,TolAbs,TolD

    TYPE(RunOptionsVar)                     :: RunOptions
    CHARACTER(LEN=200)                      :: EventName,ClassName
    PetscBool                               :: Verbose
    PetscLogEvent                           :: Event
    PetscClassId                            :: ClassID
    PetscLogDouble                          :: EventFlops=0.D0


    ClassName="System"
    CALL PetscClassIdRegister(ClassName,ClassID,ierr)
    EventName="SolveSystem"
    CALL PetscLogEventRegister(EventName,ClassID,Event,ierr)
    CALL PetscLogEventBegin(Event,PETSC_NULL_OBJECT,                 &
      & PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,      &
      & "["//ADJUSTL(TRIM(EventName))//" Event] Inizialited\n",ierr)
    

    CALL KSPCreate(PETSC_COMM_WORLD,Solver,ierr)

    CALL DMCreateGlobalVector(Gmtry%DataMngr,b,ierr)
    CALL GetInitSol(Gmtry,x,ierr)

    CALL MatSetOption(A,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr)
    CALL GetRunOptions(RunOptions,ierr)
    
    If (RunOptions%Time) THEN ! Transitory

        IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,  &
          & "["//ADJUSTL(TRIM(EventName))//                          &
          & " Event] Saving initial solution finalized\n",ierr)
        ViewName="Sol_0";ViewName=ADJUSTL(ViewName)
        CALL ViewSolution(x,ViewName,EventName,ierr)

!         CALL VecDuplicate(x,diagA,ierr)
!         CALL MatGetDiagonal(A,diagA,ierr)
        CALL MatDuplicate(A,MAT_COPY_VALUES,copyA,ierr)

        Count=1
        DO i=1,BCFld%SizeTimeZone
            DO j=1,BCFld%TimeZone(i)%SizeTime
                
                WRITE(Chari,*)i
                WRITE(Charj,*)j

                CharMsg="["//TRIM(ADJUSTL(Chari))//":"//             &
                  & TRIM(ADJUSTL(Charj))//"]"
                IF (Verbose) THEN
                  CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,     &
                    & "["//ADJUSTL(TRIM(EventName))//                &
                    & " Event] Transitory iteration "//              &
                    & TRIM(CharMsg)// " inizialited\n",ierr)
                END IF

                ! UpdateTplgy en la práctica no sirve para nada, pues 
                ! los BC se aplican con los IS del BC, no con la 
                ! geometria aún así, puede ser útil para visualizar 
                ! el cambio de la geometría en el tiempo.
!                 CALL UpdateTplgy(Gmtry,BCFld%DirichIS(i),            &
!                   & BCFld%SourceIS(i),BCFld%CauchyIS(i),ierr)

                CALL VecZeroEntries(b,ierr)
                CALL ApplyTimeDiff(PptFld,BCFld,i,j,A,b,x,ierr)
                CALL ApplyDirichlet(BCFld,i,A,b,ierr)
                CALL ApplySource(BCFld,i,b,ierr)
                CALL ApplyCauchy(BCFld,i,A,b,ierr)

                CALL KSPSetOperators(Solver,A,A,ierr)
                CALL KSPSetTolerances(Solver,PETSC_DEFAULT_REAL,     &
                  & PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,           &
                  & PETSC_DEFAULT_INTEGER,ierr)
                CALL KSPSetFromOptions(Solver,ierr)
                
                CALL KSPGetTolerances(Solver,TolR,TolAbs,TolD,MaxIts,&
                  & ierr)
                IF (Verbose)                                         &
                  & CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,   &
                  & "["//ADJUSTL(TRIM(EventName))//                  &
                  & " Event] PETSc solver monitor:\n",ierr)
                WRITE(CharTolR,*)TolR
                IF (Verbose)                                         &
                  & CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,   &
                  & "["//ADJUSTL(TRIM(EventName))//                  &
                  & " Event] Relative convergence tolerance: "//     &
                  & TRIM(CharTolR)//"\n",ierr)
                WRITE(CharTolAbs,*)TolAbs
                IF (Verbose)                                         &
                  & CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,   &
                  & "["//ADJUSTL(TRIM(EventName))//                  &
                  & " Event] Absolute convergence tolerance: "//     &
                  & TRIM(CharTolAbs)//"\n",ierr)
                WRITE(CharTolD,*)TolD
                IF (Verbose)                                         &
                  & CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,   &
                  & "["//ADJUSTL(TRIM(EventName))//                  &
                  & " Event] Divergence tolerance: "//               &
                  & TRIM(CharTolD)//"\n",ierr)
                WRITE(CharMaxIts,*)MaxIts
                IF (Verbose)                                         &
                  & CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,   &
                  & "["//ADJUSTL(TRIM(EventName))//                  &
                  & " Event] Maximun number of iterations: "//       &
                  & TRIM(CharMaxIts)//"\n",ierr)

                CALL KSPGetType(Solver,SolverType,ierr)
                IF (Verbose)                                         &
                  & CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,   &
                  & "["//ADJUSTL(TRIM(EventName))//                  &
                  & " Event] Krylov method: "//TRIM(SolverType)//    &
                  & "\n",ierr)

                CALL KSPGetPC(Solver,SolverPC,ierr)
                CALL PCGetType(SolverPC,SolverPCType,ierr)
                IF (Verbose)                                         &
                  & CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,   &
                  & "["//ADJUSTL(TRIM(EventName))//                  &
                  & " Event] Precondition method: "//                &
                  & TRIM(SolverPCType)//"\n",ierr)
                
                CALL KSPSetUp(Solver,ierr)
                CALL KSPSolve(Solver,b,x,ierr)

                CALL KSPGetConvergedReason(Solver,                   &
                  & SolverConvergedReason,ierr)
                IF (SolverConvergedReason.EQ.1) THEN
                   CharSolverConvergedReason="KSP_CONVERGED_RTOL_N"//&
                   & "ORMAL"
                ELSEIF (SolverConvergedReason.EQ.9) THEN
                   CharSolverConvergedReason="KSP_CONVERGED_ATOL_N"//&
                   & "ORMAL"
                ELSEIF (SolverConvergedReason.EQ.2) THEN
                   CharSolverConvergedReason="KSP_CONVERGED_RTOL"
                ELSEIF (SolverConvergedReason.EQ.3) THEN
                   CharSolverConvergedReason="KSP_CONVERGED_ATOL"
                ELSEIF (SolverConvergedReason.EQ.4) THEN
                   CharSolverConvergedReason="KSP_CONVERGED_ITS"
                ELSEIF (SolverConvergedReason.EQ.5) THEN
                   CharSolverConvergedReason="KSP_CONVERGED_CG_NEG"//&
                   & "_CURVE"
                ELSEIF (SolverConvergedReason.EQ.6) THEN
                   CharSolverConvergedReason="KSP_CONVERGED_CG_CON"//&
                   & "STRAINED"
                ELSEIF (SolverConvergedReason.EQ.7) THEN
                   CharSolverConvergedReason="KSP_CONVERGED_STEP_L"//&
                   & "ENGTH"
                ELSEIF (SolverConvergedReason.EQ.8) THEN
                   CharSolverConvergedReason="KSP_CONVERGED_HAPPY_"//&
                   & "BREAKDOWN"
                ELSEIF (SolverConvergedReason.EQ.-2) THEN
                   CharSolverConvergedReason="KSP_DIVERGED_NULL"
                ELSEIF (SolverConvergedReason.EQ.-3) THEN
                   CharSolverConvergedReason="KSP_DIVERGED_ITS"
                ELSEIF (SolverConvergedReason.EQ.-4) THEN
                   CharSolverConvergedReason="KSP_DIVERGED_DTOL"
                ELSEIF (SolverConvergedReason.EQ.-5) THEN
                   CharSolverConvergedReason="KSP_DIVERGED_BREAKDO"//&
                   & "WN"
                ELSEIF (SolverConvergedReason.EQ.-6) THEN
                   CharSolverConvergedReason="KSP_DIVERGED_BREAKDO"//&
                   & "WN_BICG"
                ELSEIF (SolverConvergedReason.EQ.-7) THEN
                   CharSolverConvergedReason="KSP_DIVERGED_NONSYMM"//&
                   & "ETRIC"
                ELSEIF (SolverConvergedReason.EQ.-8) THEN
                   CharSolverConvergedReason="KSP_DIVERGED_INDEFIN"//&
                   & "ITE_PC"
                ELSEIF (SolverConvergedReason.EQ.-9) THEN
                   CharSolverConvergedReason="KSP_DIVERGED_NANORINF"
                ELSEIF (SolverConvergedReason.EQ.-10) THEN
                   CharSolverConvergedReason="KSP_DIVERGED_INDEFIN"//&
                   & "ITE_MAT"
                ELSEIF (SolverConvergedReason.EQ.-11) THEN
                   CharSolverConvergedReason="KSP_DIVERGED_PCSETUP"//&
                   & "_FAILED"
                ELSEIF (SolverConvergedReason.EQ.0) THEN
                   CharSolverConvergedReason="KSP_CONVERGED_ITERATING"
                ELSE
                   CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,    &
                    & "ERROR: In Converged Reason of KSP",ierr)
                   STOP
                END IF
                
                IF (Verbose)                                         &
                  & CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,   &
                  & "["//ADJUSTL(TRIM(EventName))//                  &
                  & " Event] Converged reason: "//                   &
                  & TRIM(CharSolverconvergedreason)//"\n",ierr)
                
                IF (SolverConvergedReason.LT.0) THEN
                  CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,     &
                  & "["//ADJUSTL(TRIM(EventName))//                  &
                  & " Event] ERROR: Linear solver hasn't converged"//&
                  ".\n",ierr)
                  STOP
                END IF

                WRITE(CharCount,*)Count
                ViewName="Sol_"//TRIM(ADJUSTL(CharCount))
                ViewName=ADJUSTL(ViewName)
                CALL ViewSolution(x,ViewName,EventName,ierr)
                Count=Count+1

!                 CALL MatDiagonalSet(A,diagA,INSERT_VALUES,ierr)
                CALL MatCopy(copyA,A,SAME_NONZERO_PATTERN,ierr)

                IF (Verbose)                                         &
                  & CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,   &
                  & "["//ADJUSTL(TRIM(EventName))//                  &
                  & " Event] Transitory iteration "//TRIM(CharMsg)// &
                  & " finalized\n",ierr)

                ! Re-create the solver just to close any viewer opened 
                ! before.
                CALL KSPDestroy(Solver,ierr)
                CALL KSPCreate(PETSC_COMM_WORLD,Solver,ierr)

            END DO
        END DO
!         CALL VecDestroy(diagA,ierr)
        CALL MatDestroy(copyA,ierr)
    ELSE ! Steady

        IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,  &
          & "["//ADJUSTL(TRIM(EventName))//" Event] Steady solutio"//&
          & "n inizialited\n",ierr)

        ! UpdateTplgy en la práctica no sirve para nada, pues los BC 
        ! se aplican con los IS del BC, no con la geometria aún así, 
        ! puede ser útil para visualizar el cambio de la geometría en 
        ! el tiempo
!         CALL UpdateTplgy(Gmtry,BCFld%DirichIS(1),BCFld%SourceIS(1),  &
!           & BCFld%CauchyIS(1),ierr)

        CALL VecZeroEntries(b,ierr)
        CALL ApplyDirichlet(BCFld,1,A,b,ierr)
        CALL ApplySource(BCFld,1,b,ierr)
!         CALL ApplyCauchy(BCFld,1,A,b,ierr)


        CALL KSPSetOperators(Solver,A,A,ierr)
        CALL KSPSetTolerances(Solver,PETSC_DEFAULT_REAL,             &
          & PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,                   &
          & PETSC_DEFAULT_INTEGER,ierr)
        CALL KSPSetFromOptions(Solver,ierr)

        CALL KSPGetTolerances(Solver,TolR,TolAbs,TolD,MaxIts,ierr)
        IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,  &
          & "["//ADJUSTL(TRIM(EventName))//                          &
          & " Event] PETSc solver monitor:\n",ierr)
        WRITE(CharTolR,*)TolR
        IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,  &
          & "["//ADJUSTL(TRIM(EventName))//                          &
          & " Event] Relative convergence tolerance: "//             &
          & TRIM(CharTolR)//"\n",ierr)
        WRITE(CharTolAbs,*)TolAbs
        IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,  &
          & "["//ADJUSTL(TRIM(EventName))//                          &
          & " Event] Absolute convergence tolerance: "//             &
          & TRIM(CharTolAbs)//"\n",ierr)
        WRITE(CharTolD,*)TolD
        IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,  &
          & "["//ADJUSTL(TRIM(EventName))//                          &
          & " Event] Divergence tolerance: "//TRIM(CharTolD)//       &
          & "\n",ierr)
        WRITE(CharMaxIts,*)MaxIts
        IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,  &
          & "["//ADJUSTL(TRIM(EventName))//                          &
          & " Event] Maximun number of iterations: "//               &
          & TRIM(CharMaxIts)//"\n",ierr)

        CALL KSPGetType(Solver,SolverType,ierr)
        IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,  &
          & "["//ADJUSTL(TRIM(EventName))//                          &
          & " Event] Krylov method: "//TRIM(SolverType)//"\n",ierr)

        CALL KSPGetPC(Solver,SolverPC,ierr)
        CALL PCGetType(SolverPC,SolverPCType,ierr)
        IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,  &
          & "["//ADJUSTL(TRIM(EventName))//                          &
          & " Event] Precondition method: "//                        &
          & TRIM(SolverPCType)//"\n",ierr)
                
        CALL KSPSetUp(Solver,ierr)
        CALL KSPSolve(Solver,b,x,ierr)

        CALL KSPGetConvergedReason(Solver,SolverConvergedReason,ierr)
        IF (SolverConvergedReason.EQ.1) THEN
           CharSolverConvergedReason="KSP_CONVERGED_RTOL_NORMAL"
        ELSEIF (SolverConvergedReason.EQ.9) THEN
           CharSolverConvergedReason="KSP_CONVERGED_ATOL_NORMAL"
        ELSEIF (SolverConvergedReason.EQ.2) THEN
           CharSolverConvergedReason="KSP_CONVERGED_RTOL"
        ELSEIF (SolverConvergedReason.EQ.3) THEN
           CharSolverConvergedReason="KSP_CONVERGED_ATOL"
        ELSEIF (SolverConvergedReason.EQ.4) THEN
           CharSolverConvergedReason="KSP_CONVERGED_ITS"
        ELSEIF (SolverConvergedReason.EQ.5) THEN
           CharSolverConvergedReason="KSP_CONVERGED_CG_NEG_CURVE"
        ELSEIF (SolverConvergedReason.EQ.6) THEN
           CharSolverConvergedReason="KSP_CONVERGED_CG_CONSTRAINED"
        ELSEIF (SolverConvergedReason.EQ.7) THEN
           CharSolverConvergedReason="KSP_CONVERGED_STEP_LENGTH"
        ELSEIF (SolverConvergedReason.EQ.8) THEN
           CharSolverConvergedReason="KSP_CONVERGED_HAPPY_BREAKDOWN"
        ELSEIF (SolverConvergedReason.EQ.-2) THEN
           CharSolverConvergedReason="KSP_DIVERGED_NULL"
        ELSEIF (SolverConvergedReason.EQ.-3) THEN
           CharSolverConvergedReason="KSP_DIVERGED_ITS"
        ELSEIF (SolverConvergedReason.EQ.-4) THEN
           CharSolverConvergedReason="KSP_DIVERGED_DTOL"
        ELSEIF (SolverConvergedReason.EQ.-5) THEN
           CharSolverConvergedReason="KSP_DIVERGED_BREAKDOWN"
        ELSEIF (SolverConvergedReason.EQ.-6) THEN
           CharSolverConvergedReason="KSP_DIVERGED_BREAKDOWN_BICG"
        ELSEIF (SolverConvergedReason.EQ.-7) THEN
           CharSolverConvergedReason="KSP_DIVERGED_NONSYMMETRIC"
        ELSEIF (SolverConvergedReason.EQ.-8) THEN
           CharSolverConvergedReason="KSP_DIVERGED_INDEFINITE_PC"
        ELSEIF (SolverConvergedReason.EQ.-9) THEN
           CharSolverConvergedReason="KSP_DIVERGED_NANORINF"
        ELSEIF (SolverConvergedReason.EQ.-10) THEN
           CharSolverConvergedReason="KSP_DIVERGED_INDEFINITE_MAT"
        ELSEIF (SolverConvergedReason.EQ.-11) THEN
           CharSolverConvergedReason="KSP_DIVERGED_PCSETUP_FAILED"
        ELSEIF (SolverConvergedReason.EQ.0) THEN
           CharSolverConvergedReason="KSP_CONVERGED_ITERATING"
        ELSE
           CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,            &
            & "ERROR: In Concerged Reason of KSP",ierr)
           STOP
        END IF
        
        IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,  &
          & "["//ADJUSTL(TRIM(EventName))//                          &
          & " Event] Converged reason: "//                           &
          & TRIM(CharSolverconvergedreason)//"\n",ierr)
                

        ViewName="Sol"
        CALL ViewSolution(x,ViewName,EventName,ierr)

        IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,  &
          & "["//ADJUSTL(TRIM(EventName))//                          &
          & " Event] Steady solution finalized\n",ierr)

    END IF

    CALL KSPDestroy(Solver,ierr)

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,      &
      & "["//ADJUSTL(TRIM(EventName))//" Event] Finalized\n",ierr)
    
    CALL PetscLogFlops(EventFlops,ierr)
    CALL PetscLogEventEnd(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT, &
      & PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)


END SUBROUTINE SolveSystem

END MODULE ANISOFLOW_Solver
