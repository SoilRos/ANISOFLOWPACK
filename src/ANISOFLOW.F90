! mpiexec -n 2 ./ANISOFLOW.exe -Input_type 1 -Input_dir ../../../../DatosHet/ -Input_file_gmtry tsim_USMH.asc -Input_file_cvt matrix.mprops -Input_file_cvt_type tsim_USMH.asc -Input_file_bc_steady grid_400_400.nch_nprop_list.lateral_boundary
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

! #include <petsc/finclude/petscviewer.h>

    PetscErrorCode          :: ierr
    TYPE(Geometry)          :: Gmtry
    TYPE(PropertyField)     :: PptFld
    TYPE(BoundaryConditions):: BCFld
    Mat                     :: A
    Vec                     :: b,x

    ! PetscViewer             :: H5viewer

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Initialize program
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    CALL PetscInitialize(PETSC_COMM_WORLD,ierr)

    CALL GetGeometry(Gmtry,ierr)
    CALL GetProrperties(Gmtry,PptFld,ierr)
    CALL GetBC(Gmtry,BCFld,ierr)
    CALL GetSystem(Gmtry,PptFld,BCFld,1,A,b,x,ierr)
    CALL SolveSystem(Gmtry,BCFld,A,b,x,ierr)
    

    ! Create the HDF5 viewer
    ! CALL PetscViewerHDF5Open(PETSC_COMM_WORLD,"ANISOFLOW.h5",FILE_MODE_WRITE,H5viewer,ierr)
    ! ! CALL PetscViewerSetFromOptions(H5viewer,ierr)

    ! ! Write the H5 file 
    ! CALL VecView(x,H5viewer,ierr)

    ! ! Close the viewer
    ! CALL PetscViewerDestroy(H5viewer,ierr)

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