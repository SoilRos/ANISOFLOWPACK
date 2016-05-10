MODULE ANISOFLOW_Solver

    IMPLICIT NONE

CONTAINS

SUBROUTINE SolveSystem(Gmtry,PptFld,BCFld,A,b,x,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,BoundaryConditions,PropertyField
    USE ANISOFLOW_BuildSystem, ONLY : GetSystem
    USE ANISOFLOW_Interface

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

    PetscInt                                :: i
    PetscViewer             :: H5viewer
    CHARACTER(LEN=200)                  :: Name,tx
    PetscReal                           :: zero=0.0

    CALL KSPCreate(PETSC_COMM_WORLD,Solver,ierr)
    CALL GetRunOptions(RunOptions,ierr)

    IF (RunOptions%Time) THEN

        CALL VecDuplicate(x,diagA,ierr)
        CALL MatGetDiagonal(A,diagA,ierr)
        DO i=2,500                                               ! TEST!!!!
            PRINT*,i
            CALL GetSystem(Gmtry,PptFld,BCFld,1,i,A,b,x,ierr)   ! TEST!!!!
            CALL KSPSetOperators(Solver,A,A,ierr)
            CALL KSPSetTolerances(Solver,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,    &
                & PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
            CALL KSPSetFromOptions(Solver,ierr)
            CALL KSPSolve(Solver,b,x,ierr)

            CALL MatDiagonalSet(A,diagA,INSERT_VALUES,ierr)
            CALL VecSet(b,zero,ierr)
    ! Create the HDF5 viewer
        WRITE(tx,'(I5)')i-1
        tx=ADJUSTL(tx); tx=TRIM(tx)
    Name="ANISOFLOW_x"//TRIM(tx)//".h5"
    CALL PetscObjectSetName(x,"Altura piezometrica",ierr)
    CALL PetscViewerHDF5Open(PETSC_COMM_WORLD,Name,FILE_MODE_WRITE,H5viewer,ierr)
    ! CALL PetscViewerSetFromOptions(H5viewer,ierr)

    ! Write the H5 file 
    CALL VecView(x,H5viewer,ierr)

    ! Close the viewer
    CALL PetscViewerDestroy(H5viewer,ierr)

        END DO                                                  ! TEST!!!!

        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,             &
            & "ERROR: Transitory is not implemented yet, please use steady simulation for now",ierr)
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
