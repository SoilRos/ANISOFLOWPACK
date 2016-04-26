! ANISOFLOW_Types is a model which contains the basic data structures to the program.

MODULE ANISOFLOW_Types
    IMPLICIT NONE
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscis.h>
#include <petsc/finclude/petscdm.h>

 !  - Gemoetry: Is a data structure which manage every information needed related to geometry.
 !      + DstrMngr: Is a PETSc manager data for structured grid in 3 dimensions, provide general 
 !                  information of the grid to be managed in parallel programing.
 !      + Tplgy: Is a PETSc vector type produced by DstrMngr which contains an identificator of the 
 !               type of topology on each cell.
 !          0: Inactive cell.
 !          1: Active cell.
 !          2: Dirichlet boundary contition cell.
 !          3: Neumman on x boundary condition cell (dh/dx=0).
 !          4: Neumman on y boundary condition cell (dh/dy=0).
 !          5: Neumman on z boundary condition cell (dh/dz=0).
 !      + DirichIS: Is a PETSc index set which have a map between Dirichlet information and vecs 
 !                  produced by DstrMngr.
 !      + CauchyIS: Is a PETSc index set which have a map between Cauchy information and vecs 
 !                  produced by DstrMngr.
 !      + NeummanIS: Is a PETSc index set which have a map between Neumman information and vecs 
 !                   produced by DstrMngr.
 !      * NOTES: The variables DirichIS, CauchyIS and NeummanIS have redundant information which is 
 !             contained in the Tplgy but anyway is needed to transfer the information.

    TYPE Geometry
        DM                              :: DstrMngr
        Vec                             :: Tplgy
        IS                              :: DirichIS,CauchyIS,NeummanIS
    END TYPE Geometry


 !  - Tensor: Is a data structure which represent a full tensor, nine components, of any property.
 !      * NOTES: To use this structure is first needed fill the xx, yy, zz, xy, xz, and yz 
 !               components and then use TargetFullTensor subroutine. This process is required over 
 !               every tensor created; thereafter, is guaranteed the symmetry on the tensor.

    TYPE Tensor
        PetscReal                       :: xx,yy,zz,xy,xz,yz
        PetscReal,POINTER               :: yx,zx,zy
    END TYPE Tensor

 !  - Position: Is a data structure which describe a spatial position.
 !      + i,j,k: Integer which describe a global position on x, y, and z axes respectively.

    TYPE Position
        PetscInt                        :: i,j,k
    END TYPE Position

 !  - Property: Is a data structure wich describe completly a cell on the geometry.
 !      + Pstn: Is a Position structure which describe the global position on the grid.
 !      + StencilType: Is an integer used to describe the type of the stencil used on the model.
 !          0: Not defined.
 !          1: Star stencil. Stencil based on 6 neighbors.
 !          2: Partial box stencil. Stencil based on 18 neighbors. 
 !          3: Box stencil. Stencil based on 26 neighbors.
 !      + StencilTplgy: Is an Array of integer which contains an identificator of the type of the
 !                      topology on each cell of the stancil. The array is ordered numbering the 
 !                      cells of the stencil from upper to lowest layer, then from the lowest to the
 !                      highest values on the x and on the y axis respectively.
 !          0: Inactive cell.
 !          1: Active cell.
 !          2: Dirichlet boundary contition cell.
 !          3: Neumman on x boundary condition cell (dh/dx=0).
 !          4: Neumman on y boundary condition cell (dh/dy=0).
 !          5: Neumman on z boundary condition cell (dh/dz=0).
 !      + dx,dy,dz: A real which describe the size of the central cell in the directions on x, y, and
 !                  z axis.
 !      + dxB,dyB,dzB: A real which describe the size of the backward cell in the directions on x, 
 !                  y, and z axis.
 !      + dxF,dyF,dzF: A real which describe the size of the forward cell in the directions on x, 
 !                  y, and z axis.
 !      + CvtInterface: 
 !      + CvtCell:
 !      + CvtBx,CvtBy,CvtBz:
 !      + CvtFx,CvtFy,CvtFz:
 !      * NOTES: Do not use this structure to define each cell on a field property, it is because
 !               this structure has a redundant data that another one already has too. Instead, 
 !               use this one as temporal structure to have everything as you need of the cell in 
 !               hand.

    TYPE Property
        TYPE(Position)                  :: Pstn
        PetscInt                        :: StnclType=0
        PetscInt,ALLOCATABLE            :: StnclTplgy(:)
        PetscReal                       :: dx,dy,dz,dxB,dxF,dyB,dyF,dzB,dzF
        PetscBool                       :: CvtInterface=.TRUE.
        TYPE(Tensor)                    :: CvtCell,CvtBx,CvtFx,CvtBy,CvtFy,CvtBz,CvtFz
    END TYPE Property


    TYPE ConductivityField
        ! Conductivity field defined by few types of conductivity
        TYPE(Tensor),ALLOCATABLE                :: CvtArray(:)
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

    TYPE RunOptionsVar
        PetscBool                       :: Time=.FALSE.
        PetscInt                        :: Scheme=1
    END TYPE RunOptionsVar

    TYPE InputTypeVar
        PetscInt                        :: Gmtry,Tplgy,Cvt,SteadyBC,TransientBC
    END TYPE InputTypeVar

CONTAINS

SUBROUTINE TargetFullTensor(Tens)

    IMPLICIT NONE

    TYPE(Tensor),INTENT(INOUT)  :: Tens

    Tens%yx => TargMirrorValue(Tens%xy)
    Tens%zx => TargMirrorValue(Tens%xz)
    Tens%zy => TargMirrorValue(Tens%yz)

END SUBROUTINE TargetFullTensor

FUNCTION TargMirrorValue(Value)

    IMPLICIT NONE

    PetscReal,INTENT(IN),TARGET     :: Value
    PetscReal,POINTER               :: TargMirrorValue
    TargMirrorValue => Value

end FUNCTION TargMirrorValue

END MODULE ANISOFLOW_Types