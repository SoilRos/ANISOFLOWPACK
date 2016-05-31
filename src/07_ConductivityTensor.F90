MODULE ANISOFLOW_SetUpSystem

    IMPLICIT NONE

CONTAINS

!Subrutina ppal modulo

SUBROUTINE GetSetUpSystem(Gmtry,PptFld,BCFld,ierr)
    USE ANISOFLOW_Types
    USE ANISOFLOW_Geometry

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)               :: ierr
    TYPE(Geometry),INTENT(INOUT)               :: Gmtry

    TYPE(PropertyField),INTENT(INOUT)          :: PptFld
    TYPE(BoundaryConditions),INTENT(INOUT)     :: BCFld

    TYPE(RunOptionsVar)                        :: RunOptions
    TYPE(Geometry)                             :: UpsclGmtry

    ! It obtains run options.
    CALL GetRunOptions(RunOptions,ierr)

    IF (RunOptions%Scheme.EQ.3) THEN
        CALL GetGeometry(UpsclGmtry,ierr)
        ! Operaciones con el UpsclGmtry

        Gmtry=UpsclGmtry

    END IF


END SUBROUTINE GetSetUpSystem




END MODULE ANISOFLOW_SetUpSystem

