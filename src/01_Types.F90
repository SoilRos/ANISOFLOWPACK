
MODULE ANISOFLOW_Types
    IMPLICIT NONE
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscis.h>
#include <petsc/finclude/petscDM.h>

    TYPE Geometry
        DM                              :: DstrMngr
        Vec                             :: Tplgy
        IS                              :: DirichIS,CauchyIS,NeummanIS
    END TYPE Geometry

    TYPE RunOptionsVar
        PetscBool                       :: Time=.FALSE.
        PetscInt                        :: Scheme=1
    END TYPE RunOptionsVar

    TYPE InputTypeVar
        PetscInt                        :: Gmtry,Tplgy,Cvt,SteadyBC,TransientBC
    END TYPE InputTypeVar

    TYPE ConductivityTensor
        PetscReal                       :: xx,yy,zz,xy,xz,yz
        PetscReal,POINTER               :: yx,zx,zy
    END TYPE ConductivityTensor


    TYPE Position
        PetscInt                        :: i,j,k
    END TYPE Position

    TYPE Property
        TYPE(ConductivityTensor)        :: CvtBx,CvtFx,CvtBy,CvtFy,CvtBz,CvtFz    ! Direction (x,y,z) followed by Forward or Backward to determine position on the cube.
        PetscInt,ALLOCATABLE            :: Tplgy(:)
        TYPE(Position)                  :: Pstn
        PetscReal                       :: dx,dy,dz,dxB,dxF,dyB,dyF,dzB,dzF
        PetscBool                       :: Interface=.TRUE.     ! Every stencil must use interface conductivity
    END TYPE Property

    TYPE ConductivityField
        ! Conductivity field defined by few types of conductivity
        TYPE(ConductivityTensor),ALLOCATABLE    :: CvtArray(:)
        Vec                                     :: CvtType
        ! Complete conductivity field defined
        Vec                                     :: xxVec,yyVec,zzVec
        ! Conductivity can be saved as block or interface value, but property local must be saved as interface value
        PetscBool                               :: Interface=.FALSE.
    END TYPE ConductivityField

    TYPE PropertyField
        TYPE(ConductivityField)         :: Cvt
    END TYPE PropertyField

    TYPE StepBC
        Vec                             :: Dirich,Cauchy
    END TYPE StepBC

    TYPE BoundaryConditions
        TYPE(StepBC),ALLOCATABLE        :: Step(:)
    END TYPE BoundaryConditions

    TYPE StencilVar
        MatStencil,ALLOCATABLE          :: idx_rws(:,:),idx_clmns(:,:)
        PetscReal,ALLOCATABLE           :: Values(:)
        PetscInt                        :: Size
    END TYPE StencilVar

CONTAINS

END MODULE ANISOFLOW_Types