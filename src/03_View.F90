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
        CALL ViewVecProperty_1(x,Name,ierr)
    ELSE IF (OuputType%Sol.EQ.2) THEN
        CALL ViewVecProperty_2(x,Name,ierr)
    ELSE IF (OuputType%Sol.EQ.3) THEN
         CALL ViewVecProperty_3(x,Name,ierr)
    END IF

END SUBROUTINE ViewSolution

SUBROUTINE ViewTopology(x,Name,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    CHARACTER(LEN=200),INTENT(IN)           :: Name
    Vec,INTENT(IN)                          :: x

    TYPE(OuputTypeVar)                      :: OuputType

    CALL GetOuputType(OuputType,ierr)

    IF (OuputType%Tplgy.EQ.1) THEN
        CALL ViewVecProperty_1(x,Name,ierr)
    ELSE IF (OuputType%Tplgy.EQ.2) THEN
        CALL ViewVecProperty_2(x,Name,ierr)
    ELSE IF (OuputType%Tplgy.EQ.3) THEN
         CALL ViewVecProperty_3(x,Name,ierr)
    END IF

END SUBROUTINE ViewTopology

SUBROUTINE ViewConductivity(x,Name,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    CHARACTER(LEN=200),INTENT(IN)           :: Name
    Vec,INTENT(IN)                          :: x

    TYPE(OuputTypeVar)                      :: OuputType

    CALL GetOuputType(OuputType,ierr)

    IF (OuputType%Cvt.EQ.1) THEN
        CALL ViewVecProperty_1(x,Name,ierr)
    ELSE IF (OuputType%Cvt.EQ.2) THEN
        CALL ViewVecProperty_2(x,Name,ierr)
    ELSE IF (OuputType%Cvt.EQ.3) THEN
         CALL ViewVecProperty_3(x,Name,ierr)
    END IF

END SUBROUTINE ViewConductivity

SUBROUTINE ViewVecProperty_1(x,Name,ierr)

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
    CALL PetscViewerBinaryOpen(PETSC_COMM_WORLD,Route,FILE_MODE_WRITE,Viewer,ierr)
    CALL VecView(x,Viewer,ierr)
    CALL PetscViewerDestroy(Viewer,ierr)

    CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
        & "[] Solution "//ADJUSTL(TRIM(Name)//TRIM(Ext))//" was satisfactorily saved\n",ierr)

END SUBROUTINE ViewVecProperty_1

SUBROUTINE ViewVecProperty_2(x,Name,ierr)

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
    CALL PetscViewerASCIIOpen(PETSC_COMM_WORLD,Route,Viewer,ierr)
    CALL VecView(x,Viewer,ierr)
    CALL PetscViewerDestroy(Viewer,ierr)

    CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
        & "[] Solution "//ADJUSTL(TRIM(Name)//TRIM(Ext))//" was satisfactorily saved\n",ierr)

END SUBROUTINE ViewVecProperty_2

SUBROUTINE ViewVecProperty_3(x,Name,ierr)

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
    CALL PetscViewerHDF5Open(PETSC_COMM_WORLD,Route,FILE_MODE_WRITE,Viewer,ierr)
    CALL VecView(x,Viewer,ierr)
    CALL PetscViewerDestroy(Viewer,ierr)

    CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
        & "[] Solution "//ADJUSTL(TRIM(Name)//TRIM(Ext))//" was satisfactorily saved\n",ierr)

END SUBROUTINE ViewVecProperty_3

END MODULE ANISOFLOW_View