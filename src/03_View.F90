MODULE ANISOFLOW_View
    
    USE ANISOFLOW_Types, ONLY : OutputTypeVar

    IMPLICIT NONE

CONTAINS

SUBROUTINE ViewSolution(x,Name,StageName,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    CHARACTER(LEN=200),INTENT(IN)           :: Name,StageName
    Vec,INTENT(IN)                          :: x

    TYPE(OutputTypeVar)                      :: OutputType

    CALL GetOutputType(OutputType,ierr)

    CALL PetscObjectSetName(x,"Solution",ierr)
    IF (OutputType%Sol.EQ.1) THEN
        CALL ViewVecProperty_1(x,Name,StageName,ierr)
    ELSE IF (OutputType%Sol.EQ.2) THEN
        CALL ViewVecProperty_2(x,Name,StageName,ierr)
    ELSE IF (OutputType%Sol.EQ.3) THEN
        CALL ViewVecProperty_3(x,Name,StageName,ierr)
    END IF

END SUBROUTINE ViewSolution

SUBROUTINE ViewTopology(x,Name,StageName,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    CHARACTER(LEN=200),INTENT(IN)           :: Name,StageName
    Vec,INTENT(IN)                          :: x

    TYPE(OutputTypeVar)                      :: OutputType

    CALL GetOutputType(OutputType,ierr)

    CALL PetscObjectSetName(x,"Topology",ierr)
    IF (OutputType%Tplgy.EQ.1) THEN
        CALL ViewVecProperty_1(x,Name,StageName,ierr)
    ELSE IF (OutputType%Tplgy.EQ.2) THEN
        CALL ViewVecProperty_2(x,Name,StageName,ierr)
    ELSE IF (OutputType%Tplgy.EQ.3) THEN
        CALL ViewVecProperty_3(x,Name,StageName,ierr)
    END IF

END SUBROUTINE ViewTopology

SUBROUTINE ViewConductivity(x,Name,StageName,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    CHARACTER(LEN=200),INTENT(IN)           :: Name,StageName
    Vec,INTENT(IN)                          :: x

    TYPE(OutputTypeVar)                      :: OutputType

    CALL GetOutputType(OutputType,ierr)

    CALL PetscObjectSetName(x,"Conductivity",ierr)
    IF (OutputType%Cvt.EQ.1) THEN
        CALL ViewVecProperty_1(x,Name,StageName,ierr)
    ELSE IF (OutputType%Cvt.EQ.2) THEN
        CALL ViewVecProperty_2(x,Name,StageName,ierr)
    ELSE IF (OutputType%Cvt.EQ.3) THEN
        CALL ViewVecProperty_3(x,Name,StageName,ierr)
    END IF

END SUBROUTINE ViewConductivity

SUBROUTINE ViewStorage(x,Name,StageName,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    CHARACTER(LEN=200),INTENT(IN)           :: Name,StageName
    Vec,INTENT(IN)                          :: x

    TYPE(OutputTypeVar)                      :: OutputType

    CALL GetOutputType(OutputType,ierr)

    CALL PetscObjectSetName(x,"SpecificStorage",ierr)
    IF (OutputType%Sto.EQ.1) THEN
        CALL ViewVecProperty_1(x,Name,StageName,ierr)
    ELSE IF (OutputType%Sto.EQ.2) THEN
        CALL ViewVecProperty_2(x,Name,StageName,ierr)
    ELSE IF (OutputType%Sto.EQ.3) THEN
        CALL ViewVecProperty_3(x,Name,StageName,ierr)
    END IF

END SUBROUTINE ViewStorage

SUBROUTINE ViewProperty(x,Name,StageName,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    CHARACTER(LEN=200),INTENT(IN)           :: Name,StageName
    Vec,INTENT(IN)                          :: x

    TYPE(OutputTypeVar)                      :: OutputType

    CALL GetOutputType(OutputType,ierr)

    CALL PetscObjectSetName(x,"Property",ierr)
    IF (OutputType%Ppt.EQ.1) THEN
        CALL ViewVecProperty_1(x,Name,StageName,ierr)
    ELSE IF (OutputType%Ppt.EQ.2) THEN
        CALL ViewVecProperty_2(x,Name,StageName,ierr)
    ELSE IF (OutputType%Ppt.EQ.3) THEN
        CALL ViewVecProperty_3(x,Name,StageName,ierr)
    END IF

END SUBROUTINE ViewProperty

SUBROUTINE ViewVecProperty_1(x,Name,StageName,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscviewer.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    CHARACTER(LEN=200),INTENT(IN)           :: Name,StageName

    CHARACTER(LEN=200)                      :: OutputDir,Route
    CHARACTER(LEN=4)                        :: Ext=".bin"
    Vec,INTENT(IN)                          :: x

    PetscViewer                             :: Viewer

    CALL GetOutputDir(OutputDir,ierr)
    Route=ADJUSTL(TRIM(OutputDir)//TRIM(Name)//TRIM(Ext))
    CALL PetscViewerBinaryOpen(PETSC_COMM_WORLD,Route,FILE_MODE_WRITE,Viewer,ierr)
    CALL VecView(x,Viewer,ierr)
    CALL PetscViewerDestroy(Viewer,ierr)

    CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
        & "["//ADJUSTL(TRIM(StageName))//" Event] "//ADJUSTL(TRIM(Name)//TRIM(Ext))//" was satisfactorily saved\n",ierr)

END SUBROUTINE ViewVecProperty_1

SUBROUTINE ViewVecProperty_2(x,Name,StageName,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscviewer.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    CHARACTER(LEN=200),INTENT(IN)           :: Name,StageName

    CHARACTER(LEN=200)                      :: OutputDir,Route
    CHARACTER(LEN=4)                        :: Ext=".txt"
    Vec,INTENT(IN)                          :: x

    PetscViewer                             :: Viewer

    CALL GetOutputDir(OutputDir,ierr)
    Route=ADJUSTL(TRIM(OutputDir)//TRIM(Name)//TRIM(Ext))
    CALL PetscViewerASCIIOpen(PETSC_COMM_WORLD,Route,Viewer,ierr)
    CALL VecView(x,Viewer,ierr)
    CALL PetscViewerDestroy(Viewer,ierr)

    CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
        & "["//ADJUSTL(TRIM(StageName))//" Event] "//ADJUSTL(TRIM(Name)//TRIM(Ext))//" was satisfactorily saved\n",ierr)

END SUBROUTINE ViewVecProperty_2

SUBROUTINE ViewVecProperty_3(x,Name,StageName,ierr)

    USE ANISOFLOW_Interface

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscviewer.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    CHARACTER(LEN=200),INTENT(IN)           :: Name,StageName

    CHARACTER(LEN=200)                      :: OutputDir,Route
    CHARACTER(LEN=3)                        :: Ext=".h5"
    Vec,INTENT(IN)                          :: x

    PetscViewer                             :: Viewer

#if defined(PETSC_HAVE_HDF5)
    CALL GetOutputDir(OutputDir,ierr)
    Route=ADJUSTL(TRIM(OutputDir)//TRIM(Name)//TRIM(Ext))
    CALL PetscViewerHDF5Open(PETSC_COMM_WORLD,Route,FILE_MODE_WRITE,Viewer,ierr)
    CALL VecView(x,Viewer,ierr)
    CALL PetscViewerDestroy(Viewer,ierr)

    CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
        & "["//ADJUSTL(TRIM(StageName))//" Event] "//ADJUSTL(TRIM(Name)//TRIM(Ext))//" was satisfactorily saved\n",ierr)
#else
    CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
        & "["//ADJUSTL(TRIM(StageName))//" Event] "//ADJUSTL(TRIM(Name)//TRIM(Ext))//" was not able to save. please install hdf5 with PETSc\n",ierr)
#endif
END SUBROUTINE ViewVecProperty_3

END MODULE ANISOFLOW_View