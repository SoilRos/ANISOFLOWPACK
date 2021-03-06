PROGRAM ANISOFLOW

    USE ANISOFLOW_Types,                ONLY : Geometry,             &
                                             & PropertiesField,      &
                                             & BoundaryConditions
    USE ANISOFLOW_Interface,            ONLY : GetVerbose
    USE ANISOFLOW_Geometry,             ONLY : GetGeometry,          &
                                             & DestroyGeometry
    USE ANISOFLOW_Properties,           ONLY : GetProrperties,       &
                                             & DestroyProperties
    USE ANISOFLOW_BoundaryConditions,   ONLY : GetBC,DestroyBC
    USE ANISOFLOW_BuildSystem,          ONLY : BuildSystem,          &
                                              & DestroySystem
    USE ANISOFLOW_Solver,               ONLY : SolveSystem

    IMPLICIT NONE
    
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>

    PetscErrorCode              :: ierr
    PetscBool                   :: Verbose
    TYPE(Geometry)              :: Gmtry

    TYPE(BoundaryConditions)    :: BCFld
    Mat                         :: A
    Vec                         :: b,x

    PetscLogStage               :: Stage
    CHARACTER(LEN=200)          :: StageName
    TYPE(PropertiesField)       :: PptFld


    CALL PetscInitialize(PETSC_NULL_CHARACTER,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,      &
        & "[Main Stage] Inizialited\n",ierr)

    StageName="Loading Files"
    CALL PetscLogStageRegister(StageName, stage,ierr)
    CALL PetscLogStagePush(Stage,ierr)

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,      &
        & "["//ADJUSTL(TRIM(StageName))//" Stage] Inizialited\n",ierr)

    CALL GetGeometry(PETSC_COMM_WORLD,Gmtry,ierr)
    CALL GetProrperties(Gmtry,PptFld,ierr)
    CALL GetBC(Gmtry,BCFld,ierr)

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,      &
        & "["//ADJUSTL(TRIM(StageName))//" Stage] Finalized\n",ierr)
    CALL PetscLogStagePop(Stage,ierr)

    StageName="Setting up"
    CALL PetscLogStageRegister(StageName, stage,ierr)
    CALL PetscLogStagePush(Stage,ierr)

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,      &
        & "["//ADJUSTL(TRIM(StageName))//" Stage] Inizialited\n",ierr)

    CALL BuildSystem(Gmtry,PptFld,A,ierr)

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,      &
        & "["//ADJUSTL(TRIM(StageName))//" Stage] Finalized\n",ierr)
    CALL PetscLogStagePop(Stage,ierr)

    StageName="Solving"
    CALL PetscLogStageRegister(StageName, stage,ierr)
    CALL PetscLogStagePush(Stage,ierr)

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,      &
        & "["//ADJUSTL(TRIM(StageName))//" Stage] Inizialited\n",ierr)

    CALL SolveSystem(Gmtry,PptFld,BCFld,A,b,x,ierr)

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,      &
        & "["//ADJUSTL(TRIM(StageName))//" Stage] Finalized\n",ierr)
    CALL PetscLogStagePop(Stage,ierr)

    StageName="Destroying obj"
    CALL PetscLogStageRegister(StageName, stage,ierr)
    CALL PetscLogStagePush(Stage,ierr)

    CALL DestroySystem(A,b,x,ierr)
    CALL DestroyBC(BCFld,ierr)
    CALL DestroyProperties(PptFld,ierr)
    CALL DestroyGeometry(Gmtry,ierr)

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,      &
        & "["//ADJUSTL(TRIM(StageName))//" Stage] Finalized\n",ierr)
    CALL PetscLogStagePop(Stage,ierr)

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,      &
        & "[Main Stage] Finalized\n",ierr)

    CALL PetscFinalize(ierr)

END PROGRAM

! Use PetscDataTypeToMPIDataType subroutine to Send/Recive MPI values 
! when PETSc team add a Fortran wrapper function for it.
