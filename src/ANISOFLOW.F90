PROGRAM ANISOFLOW

    USE ANISOFLOW_Types,                ONLY : Geometry,PropertyField,BoundaryConditions
    USE ANISOFLOW_Interface,            ONLY : GetVerbose
    USE ANISOFLOW_Geometry,             ONLY : GetGeometry,DestroyGeometry
    USE ANISOFLOW_Properties,           ONLY : GetProrperties,DestroyProperties
    USE ANISOFLOW_BoundaryConditions,   ONLY : GetBC,DestroyBC
    USE ANISOFLOW_BuildSystem,          ONLY : GetSystem,DestroySystem
    USE ANISOFLOW_Solver,               ONLY : SolveSystem

    IMPLICIT NONE
    
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>

    PetscErrorCode          :: ierr
    PetscBool               :: Verbose
    TYPE(Geometry)          :: Gmtry
    TYPE(PropertyField)     :: PptFld
    TYPE(BoundaryConditions):: BCFld
    Mat                     :: A
    Vec                     :: b,x

    CALL PetscInitialize(PETSC_COMM_WORLD,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[Main Stage] Inizialited\n",ierr)

    CALL GetGeometry(Gmtry,ierr)
    CALL GetProrperties(Gmtry,PptFld,ierr)
    CALL GetBC(Gmtry,BCFld,ierr)
    CALL GetSystem(Gmtry,PptFld,BCFld,1,1,A,b,x,ierr)
    CALL SolveSystem(Gmtry,PptFld,BCFld,A,b,x,ierr)

    CALL DestroySystem(A,b,x,ierr)
    CALL DestroyBC(BCFld,ierr)
    CALL DestroyProperties(PptFld,ierr)
    CALL DestroyGeometry(Gmtry,ierr)

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[Main Stage] Finalized\n",ierr)

    CALL PetscFinalize(ierr)

END PROGRAM

! Use PetscDataTypeToMPIDataType subroutine to Send/Recive MPI values when PETSc team add 
! a Fortran wrapper function for it.