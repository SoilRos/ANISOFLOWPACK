MODULE ANISOFLOW_Types

! ANISOFLOW_Types it's a module that contains the basic data structures to ANISOFLOWPACK.

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscis.h>
#include <petsc/finclude/petscdm.h>

 !  - Geometry: It's a data structure that manages every information needed related to geometry.
 !    > VARIABLES: DataMngr, Tplgy, DirichIS, CauchyIS, NeummanIS.
 !      + DataMngr: It's a PETSc manager data for structured grid in 3 dimensions, providing general 
 !                  information on the grid to be managed in parallel programming.
 !      + Tplgy: It's a PETSc vector type produced by DataMngr that contains a topology identifier
 !               in each cell.
 !          0: Inactive cell.
 !          1: Active cell.
 !          2: Dirichlet boundary condition cell.
 !          3: Neumman on x boundary condition cell (dh/dx=0).
 !          4: Neumman on y boundary condition cell (dh/dy=0).
 !          5: Neumman on z boundary condition cell (dh/dz=0).
 !      + DirichIS: It's a PETSc index set that has a map between Dirichlet information and vecs 
 !                  produced by DataMngr.
 !      + CauchyIS: It's a PETSc index set that has a map between Cauchy information and vecs 
 !                  produced by DataMngr.
 !      + NeummanIS: It's a PETSc index set that has a map between Neumman information and vecs 
 !                   produced by DataMngr.
 !    > NOTES: The variables DirichIS, CauchyIS and NeummanIS have redundant information  that 
 !             contains the Tplgy. In any case is needed to transfer the information.

    TYPE Geometry
        DM                              :: DataMngr
        Vec                             :: Tplgy
        IS                              :: DirichIS,CauchyIS,NeummanIS
    END TYPE Geometry


 !  - Tensor: It's a data structure that represents a full tensor, nine components, of any property.
 !    > VARIABLES: xx, yy, zz, xy, xz, yz, yx, zx, zy.
 !      + xx,yy,zz,xy,xz,yz,yx,zx,zy: Conductivity components.
 !    > NOTES: To use this structure is first needed to fill the xx, yy, zz, xy, xz, and yz 
 !               components and then use TargetFullTensor subroutine to have (xy=yx, xz=zx, and 
 !               yz=zy). This process is required over every tensor created; after that is 
 !               guaranteed the symmetry on the tensor.

    TYPE Tensor
        PetscReal                       :: xx,yy,zz,xy,xz,yz
        PetscReal,POINTER               :: yx,zx,zy
    END TYPE Tensor

 !  - Position: It's a data structure that describes a spatial position.
 !      + i,j,k: It's an integer that describes a global position on x, y, and z axes respectively.

    TYPE Position
        PetscInt                        :: i,j,k
    END TYPE Position

 !  - Property: It's a data structure that describes completely a cell on the geometry.
 !    > VARIABLES: Pstn, StencilType, StencilTplgy, dx, dy, dz, dxB, dyB, dzB, dxF, dyF, dzF,
 !                 CvtOnBlock, CvtBlock, CvtOnInterface, CvtBx, CvtBy, CvtBz, CvtFx, CvtFy, CvtFz.
 !      + Pstn: It's a Position structure that describes the global position on the grid.
 !      + StencilType: It's an integer used to describe the type of the stencil used on the model.
 !          0: Not defined. (Default)
 !          1: Star stencil. Stencil based on 6 neighbors.
 !          2: Partial box stencil. Stencil based on 18 neighbors. 
 !          3: Box stencil. Stencil based on 26 neighbors.
 !      + StencilTplgy: It's an Array of integers that contains an identifier topology in each cell 
 !                      of the stencil. The array is ordered numbering the cells of the stencil from 
 !                      upper to lowest layer, then from the lowest to the highest values on the x 
 !                      and the y axis respectively.
 !          0: Inactive cell.
 !          1: Active cell.
 !          2: Dirichlet boundary condition cell.
 !          3: Neumman on x boundary condition cell (dh/dx=0).
 !          4: Neumman on y boundary condition cell (dh/dy=0).
 !          5: Neumman on z boundary condition cell (dh/dz=0).
 !      + dx,dy,dz: A real that describes the size of the central cell in the directions on x, y,
 !                  and z axis.
 !      + dxB,dyB,dzB: A real that describes the size of the backward cell in the directions on x, 
 !                  y, and z axis.
 !      + dxF,dyF,dzF: A real that describes the size of the forward cell in the directions on x, 
 !                  y, and z axis.
 !      + CvtOnBlock: It's a boolean that shows if the block has the conductivity represented 
 !                    on the block.
 !      + CvtBlock: It's a Tensor of conductivities of the medium of the block.
 !      + CvtOnInterface: It's a boolean that shows if the block has the conductivity represented 
 !                        on the interfaces.
 !      + CvtBx,CvtBy,CvtBz: It's a tensor of conductivities on the interface of the block in
 !                           backward direction of the cartesian axis respectively.
 !      + CvtFx,CvtFy,CvtFz: It's a tensor of conductivities on the interface of the block in  
 !                           forward direction of the cartesian axes respectively.
 !    > NOTES: Don't use this structure to define each cell on a field of properties, it is 
 !             because this structure has a redundant data that another one already has too. 
 !             Instead, use this one as a temporal structure to have everything as you need of the 
 !             cell in hand.

    TYPE Property
        TYPE(Position)                  :: Pstn
        PetscInt                        :: StnclType=0
        PetscInt,ALLOCATABLE            :: StnclTplgy(:)
        PetscReal                       :: dx,dy,dz,dxB,dxF,dyB,dyF,dzB,dzF
        PetscBool                       :: CvtOnInterface=.FALSE.,CvtOnBlock=.FALSE.
        TYPE(Tensor)                    :: CvtBlock,CvtBx,CvtFx,CvtBy,CvtFy,CvtBz,CvtFz
    END TYPE Property

 !  - ConductivityField: It's a data structure that stores the field of conductivities of every
 !                       block. This structure may store conductivities defined by field zones or
 !                       where each conductivity block has a different value.
 !    > VARIABLES: DefinedByZones, DefindeByInteface, CvtZones, CvtType, xxVec, yyVec, zzVec.
 !      + DefinedByZones: It's a boolean describing if the ConductivityField is stored by zones or 
 !                        not. If the boolean is .TRUE., it's necessary save the information in 
 !                        CvtZones and CvtYpe, otherwise it's necessary save the information in 
 !                        xxVec,yyVec,zzVec
 !      + DefindeByInteface: It's a boolean describing if the ConductivityField stores the
 !                           information on the interface or each block. If the boolean is .TRUE.,
 !                           the CvtType or xxVec, yyVec, and zzVec need an additional component
 !                           in each direction in the first processor.
 !      + CvtZones: It's an Array of Tenors that contain many Tensors of conductivity as defined 
 !                  zones. This variable needs to be stored in all processors.
 !      + CvtType: It's a PETSc vector that contains an identifier for each cell, the identifier 
 !                 value correspond to a zone value.
 !      + xxVec,yyVec,zzVec: Contain a tensor conductivity component by each cell. 

    TYPE ConductivityField
        PetscBool                               :: DefinedByZones
        PetscBool                               :: DefinedByInterface=.FALSE.
        ! Conductivity defined by zones:
        TYPE(Tensor),ALLOCATABLE                :: CvtZone(:)
        Vec                                     :: CvtType
        ! Conductivity defined on every cell:
        Vec                                     :: xxVec,yyVec,zzVec
    END TYPE ConductivityField

 !  - ProperyField: It's a data structure wich contains fields of differents properties.
 !    > VARIABLES: Cvt.
 !      + Cvt: It's a ConductivityField data structure wich contain a field of conductivity.

    TYPE PropertyField
        TYPE(ConductivityField)         :: Cvt
    END TYPE PropertyField

 !  - StepBC: It's a data structure that contains a Boundary Condition set.
 !    > VARIABLES: Diricl, Cauchy.
 !      + Diricl: It's a PETSc vector with the values of the Dirichlet cells. It must have the same 
 !                size as Dirichlet identifiers defined in Geometry data structure.
 !      + Cauchy: It's a PETSc vector with the values of the Cauchy cells. It must have the same 
 !                size as Cauchy identifiers defined in Geometry data structure.
 !    > NOTES: The Index Sets defined in Geometry for each boundary condition point the value of  
 !             Dirich and Cauchy vectors to the real position of every global vector created with
 !             DataMngr defined in Geometry.
    TYPE StepBC
        Vec                             :: Dirich,Cauchy
    END TYPE StepBC

 !  - BoundaryConditions: It's a data structure wich contains a series Boundary Conditions to be  
 !                        used in several times or zones of time.
 !      + Step: It's a StepBc array with several values of Boundary Conditions.

    TYPE BoundaryConditions
        TYPE(StepBC),ALLOCATABLE        :: Step(:)
    END TYPE BoundaryConditions

 !  - StensilVar: It's a auxilar data structure to build the matrix.
 !    > VARIABLES: idx_rws, idx_clmn, Values, Size.
 !      + idx_rws: It's a PETSc MatStencil structure that contains the information of the rows that 
 !                 will be modified in the matrix. The allocation to the first component is always 
 !                 four, to the second one it depends on how many rows are necessary to change, 
 !                 usually one.
 !      + idx_clmn: It's a PETSc MatStencil structure that contanis the information of the columns 
 !                  that will be modified on the matrix. The allocation to the first component is  
 !                  always four, to the second one it depends on how many colums are necessary to 
 !                  change.
 !      + Values: It's an array that stores the values to be added to the position indicated by 
 !                idx_rws and idx_clmns.
 !      + Size: It's an integer that says how many values will be modified on the matrix.
 !    > NOTES: An example of usage of the StencilVar of 3 values in a row is the following:
 !          
 !!          StencilVar%Size=3                           ! Amount of values of the Stencil
 !!
 !!          ALLOCATE(StencilVar%idx_rws(4,1))           ! One row to modify.
 !!          ALLOCATE(StencilVar%idx_clmn(4,3))          ! Four values of the row to modify.
 !!          ALLOCATE(StencilVar%Values(StencilVar%Size))! Amount of values of the Stencil
 !!
 !!          StencilVar%idx_rws(MatStencil_i,1) = i      ! The value of the row to be modified will 
 !!          StencilVar%idx_rws(MatStencil_j,1) = j      !   be associated with the position i, j, 
 !!          StencilVar%idx_rws(MatStencil_k,1) = k      !   and k.
 !!
 !!          StencilVar%idx_clmn(MatStencil_i,:) = i     ! The value of the columns to be modified 
 !!          StencilVar%idx_clmn(MatStencil_j,1) = j-1   !   will be associated with the position i,
 !!          StencilVar%idx_clmn(MatStencil_j,2) = j     !   j-1, and k to the fisrt (1) value, i,
 !!          StencilVar%idx_clmn(MatStencil_j,3) = j+1   !   j, and k to the second (2) value, and
 !!          StencilVar%idx_clmn(MatStencil_k,:) = k     !   i, j+1, and k to the third value.
 !!          
 !!          StencilVar%Values(1)=1.0                    ! Values of the stencil that will be added 
 !!          StencilVar%Values(2)=2.0                    !    to the positions defined with idx_rws
 !!          StencilVar%Values(3)=1.0                    !    and idx_clmn.
 !                                  
 !             Finally, the StencilVar needs to be added to the matrix with MatSetValuesStencil
 !             that is a PETSc function.

    TYPE StencilVar
        MatStencil,ALLOCATABLE          :: idx_rws(:,:),idx_clmns(:,:)
        PetscReal,ALLOCATABLE           :: Values(:)
        PetscInt                        :: Size
    END TYPE StencilVar

 !  - RunOptionsVar: It's a data structure that contains all options related to the running.
 !    > VARIABLES: Time, Scheme.
 !      + Time: It's a boolean that defines if the running is steady or transient. The default value
 !              is steady.
 !      + Scheme: It's an integer that defines the scheme stencil to solve the underground flow 
 !                equation.
 !          1: Traditional solution. MODFLOW uses it, It uses an interface conductivity determined
 !             by the harmonic. Default use if it's not indicated.
 !          2: Li model. It uses an interface conductivity to build the stencil.
 !          3: ANISOFLOWPACK model. -----

    TYPE RunOptionsVar
        PetscBool                       :: Time=.FALSE.
        PetscInt                        :: Scheme=1
    END TYPE RunOptionsVar

 !  - InputTypeVar: It's a collection of integer that defines the type of input to be used in the
 !                 program.
 !    > VARIABLES: Gmtry,Tplgy,Cvt,BC
 !      + Gmtry: It's an integer that defines a type of file that will provide the domain size and
 !               the cell center coordinates.
 !          1: Defined by Blessent. An example is provided in "../ex/Blessent/in/tsim_USMH.asc"
 !          2: Defined by Perez. An example is provided in "../ex/Perez/in/sanpck.domnRST"
 !      + Tplgy: It's an integer that defines a type of file that will provide the identifiers to
 !               each cell.
 !          1: Default topology, it doesn't need a file. Every face of the domain is a Neumman
 !             condition but on the first layer where the boundary is dirichlet.
 !      + Cvt: It's an integer that defines a type of file that will provide conductivity to each 
 !             cell or zone.
 !          1: It's a pair of files that provide block conductivities by zones and another one that
 !             provide an identifier of zone to each cell. An example is provided in "../ex/
 !             Blessent/in/matrix.mprops" and "../ex/Blessent/in/tsim_USMH.asc".
 !      + BC: It's a integer that define a type of file that will provide the boundary conditions.
 !          1: It's a file that only provide dirichlet condition and their position on the grid. An
 !             example is provided in "../Blessent/in/grid_400_400.nch_nprop_list.lateral_boundary".

    TYPE InputTypeVar
        PetscInt                        :: Gmtry,Tplgy,Cvt,BC
    END TYPE InputTypeVar

CONTAINS

 !  - TargetFullTensor: It's an auxiliary routine to guarantee the Tensor variable symetry.
 !    > IN: Tensor
 !      + Tensor: A tensor structure with xx, yy, zz, xy, xz, and yz defined.
 !    > OUT: Tensor
 !      + Tensor: A tensor with yx, zx, and zy values pointed to xy, xz, and yz values respectively.

SUBROUTINE TargetFullTensor(Tens)

    IMPLICIT NONE

    TYPE(Tensor),INTENT(INOUT)  :: Tens

    Tens%yx => TargMirrorValue(Tens%xy)
    Tens%zx => TargMirrorValue(Tens%xz)
    Tens%zy => TargMirrorValue(Tens%yz)

END SUBROUTINE TargetFullTensor

 !  - TargMirrorValue: It's a auxilar function to target a real that will be pointed later.
 !    > IN: Value
 !      + Values: It's a real.
 !    > OUT: Value
 !      + Values: It's the same input real but now as a target.

FUNCTION TargMirrorValue(Value)

    IMPLICIT NONE

    PetscReal,INTENT(IN),TARGET     :: Value
    PetscReal,POINTER               :: TargMirrorValue
    TargMirrorValue => Value

end FUNCTION TargMirrorValue

END MODULE ANISOFLOW_Types