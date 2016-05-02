PROGRAM ANISOFLOW
    USE ANISOFLOW_Geometry
    USE ANISOFLOW_Properties
    USE ANISOFLOW_BoundaryConditions
    USE ANISOFLOW_BuildSystem
    USE ANISOFLOW_Solver
    IMPLICIT NONE
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>


    PetscErrorCode          :: ierr
    TYPE(Geometry)          :: Gmtry
    TYPE(PropertyField)     :: PptFld
    TYPE(BoundaryConditions):: BCFld
    Mat                     :: A
    Vec                     :: b,x

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Initialize program
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    CALL PetscInitialize(PETSC_COMM_WORLD,ierr)

    CALL GetGeometry(Gmtry,ierr)
    CALL GetProrperties(Gmtry,PptFld,ierr)
    CALL GetBC(Gmtry,BCFld,ierr)
    CALL GetSystem(Gmtry,PptFld,BCFld,1,A,b,x,ierr)
    CALL SolveSystem(Gmtry,BCFld,A,b,x,ierr)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Finalize program
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    CALL SystemDestroy(A,b,x,ierr)
    CALL BCDestroy(BCFld,ierr)
    CALL PropertiesDestroy(PptFld,ierr)
    CALL GeometryDestroy(Gmtry,ierr)
    CALL PetscFinalize(ierr)

END PROGRAM

! Hay que preguntarle a los del PETSc con cual "derived type" se deben transmitir
! las variables para que el programa siga siendo portable. (Corregir en Properties y Geometry)