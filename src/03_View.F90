MODULE ANISOFLOW_View
    
    USE ANISOFLOW_Types, ONLY : OuputTypeVar

    IMPLICIT NONE

CONTAINS

SUBROUTINE ViewSolution(x,Name,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    CHARACTER(LEN=200),INTENT(IN)           :: Name
    Vec,INTENT(IN)                          :: x

    TYPE(OuputTypeVar)                      :: OuputType

    CALL GetOuputType(OuputType,ierr)

    IF (OuputType%Sol.EQ.1) THEN
        CALL ViewSolution_1(x,Name,ierr)
    ELSE 
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "ERROR: OuputType sol wrong\n",ierr)
        STOP
    END IF

END SUBROUTINE ViewSolution

SUBROUTINE ViewSolution_1(x,Name,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscviewer.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    CHARACTER(LEN=200),INTENT(IN)           :: Name
    CHARACTER(LEN=200)                      :: OuputDir,Route
    Vec,INTENT(IN)                          :: x

    PetscViewer                             :: H5viewer

    CALL GetOuputDir(OuputDir,ierr)
    Route=ADJUSTL(TRIM(OuputDir)//TRIM(Name))
    CALL PetscObjectSetName(x,"PiezometricHead",ierr)
    ! CALL PetscViewerHDF5Open(PETSC_COMM_WORLD,Route,FILE_MODE_WRITE,H5viewer,ierr)
    ! CALL VecView(x,H5viewer,ierr)
    ! CALL PetscViewerDestroy(H5viewer,ierr)

END SUBROUTINE ViewSolution_1

END MODULE ANISOFLOW_View