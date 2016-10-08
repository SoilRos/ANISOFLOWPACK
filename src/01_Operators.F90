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

!         TYPE(TENSOR) FUNCTION TensorArmonic(x,y)
!             USE ANISOFLOW_Types, ONLY : Tensor, TargetFullTensor
!             IMPLICIT NONE
! #include <petsc/finclude/petsc.h>
!             TYPE(Tensor), INTENT(IN) :: x,y
!         END FUNCTION TensorArmonic
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

    IF ((x==0).OR.(y==0)) THEN
        RealArmonic=0.D0
    ELSE
        RealArmonic = 2.0/(1.0/x+1.0/y)
    END IF

END FUNCTION RealArmonic


! TYPE(TENSOR) FUNCTION TensorArmonic(x,y)

!     USE ANISOFLOW_Types, ONLY : Tensor, TargetFullTensor

!     IMPLICIT NONE

! #include <petsc/finclude/petsc.h>

!     TYPE(Tensor), INTENT(IN) :: x,y

!     IF ((x%xx==0).OR.(y%xx==0)) THEN
!         TensorArmonic%xx=0.D0
!     ELSE
!         TensorArmonic%xx = 2.0/(1.0/x%xx+1.0/y%xx)
!     END IF

!     IF ((x%yy==0).OR.(y%yy==0)) THEN
!         TensorArmonic%yy=0.D0
!     ELSE
!         TensorArmonic%yy = 2.0/(1.0/x%yy+1.0/y%yy)
!     END IF

!     IF ((x%zz==0).OR.(y%zz==0)) THEN
!         TensorArmonic%zz=0.D0
!     ELSE
!         TensorArmonic%zz = 2.0/(1.0/x%zz+1.0/y%zz)
!     END IF

!     IF ((x%xy==0).OR.(y%xy==0)) THEN
!         TensorArmonic%xy=0.D0
!     ELSE
!         TensorArmonic%xy = 2.0/(1.0/x%xy+1.0/y%xy)
!     END IF

!     IF ((x%xz==0).OR.(y%xz==0)) THEN
!         TensorArmonic%xz=0.D0
!     ELSE
!         TensorArmonic%xz = 2.0/(1.0/x%xz+1.0/y%xz)
!     END IF

!     IF ((x%yz==0).OR.(y%yz==0)) THEN
!         TensorArmonic%yz=0.D0
!     ELSE
!         TensorArmonic%yz = 2.0/(1.0/x%yz+1.0/y%yz)
!     END IF

!     CALL TargetFullTensor(TensorArmonic)

! END FUNCTION TensorArmonic