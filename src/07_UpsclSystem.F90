MODULE ANISOFLOW_UpsclSystem

    IMPLICIT NONE

CONTAINS

!This subrutine manage all subrutines related with upscaling procedures
!

SUBROUTINE GetUpsclSystem(Gmtry,PptFld,BCFld,ierr)
    USE ANISOFLOW_Types
    USE ANISOFLOW_Geometry

#include <petsc/finclude/petscsys.h>


    PetscErrorCode,INTENT(INOUT)               :: ierr
    TYPE(Geometry),INTENT(INOUT)               :: Gmtry
    TYPE(PropertyField),INTENT(INOUT)          :: PptFld
    TYPE(BoundaryConditions),INTENT(INOUT)     :: BCFld

    TYPE(RunOptionsVar)                        :: RunOptions
    TYPE(Geometry)                             :: UpsclGmtry
    PetscInt                                   :: u,i,width(3),skin

    !Skin used in upscaling procedure
    PARAMETER(skin=5)

    ! It obtains run options.
    CALL GetRunOptions(RunOptions,ierr)

    IF (RunOptions%Scheme.EQ.3) THEN
        !Get geometry at coarse scale
        CALL GetGeometry(UpsclGmtry,ierr)
        ! Operaciones con el UpsclGmtry
        CALL DMDAGetInfo(UpsclGmtry%DataMngr,PETSC_NULL_INTEGER,width(1),width(2),&
        & width(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,ierr)

        DO i=1,width(1)*width(2)*width(3)


        END DO

        Gmtry=UpsclGmtry

    END IF


END SUBROUTINE GetUpsclSystem

SUBROUTINE GetUpscldBC(Gmtry,PptFld,BCFld,ierr)


END SUBROUTINE GetUpscldBC

END MODULE ANISOFLOW_SetUpsclSystem

