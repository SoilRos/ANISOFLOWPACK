MODULE ANISOFLOW_Operators

    IMPLICIT NONE

    INTERFACE ASSIGNMENT(=)
        SUBROUTINE EqualTensors(y,x)
            USE ANISOFLOW_Types, ONLY : Tensor
            IMPLICIT NONE
            TYPE(Tensor),INTENT(OUT)    :: y
            TYPE(Tensor),INTENT(IN)     :: x
        END SUBROUTINE EqualTensors
    END INTERFACE

    INTERFACE OPERATOR(.ARMONIC.)
        FUNCTION RealArmonic(x,y)
            IMPLICIT NONE
#include <petsc/finclude/petsc.h>
            PetscReal, INTENT(IN) :: x,y
            PetscReal             :: RealArmonic
        END FUNCTION RealArmonic

    END INTERFACE

END MODULE ANISOFLOW_Operators

SUBROUTINE EqualTensors(y,x)

    USE ANISOFLOW_Types, ONLY : Tensor, TargetFullTensor

    IMPLICIT NONE

    TYPE(Tensor),INTENT(OUT)    :: y
    TYPE(Tensor),INTENT(IN)     :: x

    y%xx=x%xx
    y%yy=x%yy
    y%zz=x%zz
    y%xy=x%xy
    y%xz=x%xz
    y%yz=x%yz
    CALL TargetFullTensor(y)

END SUBROUTINE EqualTensors


FUNCTION RealArmonic(x,y)

    IMPLICIT NONE

#include <petsc/finclude/petsc.h>

    PetscReal, INTENT(IN) :: x,y

    PetscReal             :: RealArmonic

!     PetscInt              :: InterfaceOperatorType
!     COMMON                /INFERFACEOPERATOR/ InterfaceOperatorType

    IF ((x==0).OR.(y==0)) THEN
        RealArmonic=0.D0
    ELSE
        RealArmonic = 2.0/(1.0/x+1.0/y)
    END IF

END FUNCTION RealArmonic
