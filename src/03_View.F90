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
    ELSE IF (OuputType%Sol.EQ.2) THEN
        CALL ViewSolution_2(x,Name,ierr)
    ! ELSE IF (OuputType%Sol.EQ.3) THEN
    !     CALL ViewSolution_3(x,Name,ierr)
    ! ELSE IF (OuputType%Sol.EQ.4) THEN
    !     CALL ViewSolution_3(x,Name,ierr)
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
    CHARACTER(LEN=4)                        :: Ext=".bin"
    Vec,INTENT(IN)                          :: x

    PetscViewer                             :: Viewer

    CALL GetOuputDir(OuputDir,ierr)
    Route=ADJUSTL(TRIM(OuputDir)//TRIM(Name)//TRIM(Ext))
    CALL PetscObjectSetName(x,"PiezometricHead",ierr)
    CALL PetscViewerBinaryOpen(PETSC_COMM_WORLD,Route,FILE_MODE_WRITE,Viewer,ierr)
    CALL VecView(x,Viewer,ierr)
    CALL PetscViewerDestroy(Viewer,ierr)

END SUBROUTINE ViewSolution_1

SUBROUTINE ViewSolution_2(x,Name,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscviewer.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    CHARACTER(LEN=200),INTENT(IN)           :: Name

    CHARACTER(LEN=200)                      :: OuputDir,Route
    CHARACTER(LEN=4)                        :: Ext=".txt"
    Vec,INTENT(IN)                          :: x

    PetscViewer                             :: Viewer

    CALL GetOuputDir(OuputDir,ierr)
    Route=ADJUSTL(TRIM(OuputDir)//TRIM(Name)//TRIM(Ext))
    CALL PetscObjectSetName(x,"PiezometricHead",ierr)
    CALL PetscViewerASCIIOpen(PETSC_COMM_WORLD,Route,Viewer,ierr)

    CALL VecView(x,Viewer,ierr)
    CALL PetscViewerDestroy(Viewer,ierr)

END SUBROUTINE ViewSolution_2

SUBROUTINE ViewSolution_3(x,Name,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscviewer.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    CHARACTER(LEN=200),INTENT(IN)           :: Name

    CHARACTER(LEN=200)                      :: OuputDir,Route
    CHARACTER(LEN=4)                        :: Ext=".mat"
    Vec,INTENT(IN)                          :: x

    PetscViewer                             :: Viewer

    CALL GetOuputDir(OuputDir,ierr)
    Route=ADJUSTL(TRIM(OuputDir)//TRIM(Name)//TRIM(Ext))
    CALL PetscObjectSetName(x,"PiezometricHead",ierr)
    ! CALL PetscViewerMatlabOpen(PETSC_COMM_WORLD,Route,FILE_MODE_WRITE,Viewer,ierr)

    CALL VecView(x,Viewer,ierr)
    CALL PetscViewerDestroy(Viewer,ierr)

END SUBROUTINE ViewSolution_3

SUBROUTINE ViewSolution_4(x,Name,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscviewer.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    CHARACTER(LEN=200),INTENT(IN)           :: Name

    CHARACTER(LEN=200)                      :: OuputDir,Route
    CHARACTER(LEN=3)                        :: Ext=".h5"
    Vec,INTENT(IN)                          :: x

    PetscViewer                             :: Viewer

    CALL GetOuputDir(OuputDir,ierr)
    Route=ADJUSTL(TRIM(OuputDir)//TRIM(Name)//TRIM(Ext))
    CALL PetscObjectSetName(x,"PiezometricHead",ierr)
    CALL PetscViewerHDF5Open(PETSC_COMM_WORLD,Route,FILE_MODE_WRITE,Viewer,ierr)
    CALL VecView(x,Viewer,ierr)
    CALL PetscViewerDestroy(Viewer,ierr)

END SUBROUTINE ViewSolution_4

END MODULE ANISOFLOW_View