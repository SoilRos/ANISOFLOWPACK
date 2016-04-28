! ANISOFLOW_Types it's a module that contains the basic data structures to the program.

MODULE ANISOFLOW_Types
    IMPLICIT NONE
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscis.h>
#include <petsc/finclude/petscdm.h>

 !  - Gemoetry: It's a data structure that manage every information needed related to geometry.
 !      + DstrMngr: It's a PETSc manager data for structured grid in 3 dimensions, provide general 
 !                  information of the grid to be managed in parallel programing.
 !      + Tplgy: It's a PETSc vector type produced by DstrMngr that contains an identificator of the 
 !               type of topology on each cell.
 !          0: Inactive cell.
 !          1: Active cell.
 !          2: Dirichlet boundary contition cell.
 !          3: Neumman on x boundary condition cell (dh/dx=0).
 !          4: Neumman on y boundary condition cell (dh/dy=0).
 !          5: Neumman on z boundary condition cell (dh/dz=0).
 !      + DirichIS: It's a PETSc index set that have a map between Dirichlet information and vecs 
 !                  produced by DstrMngr.
 !      + CauchyIS: It's a PETSc index set that have a map between Cauchy information and vecs 
 !                  produced by DstrMngr.
 !      + NeummanIS: It's a PETSc index set that have a map between Neumman information and vecs 
 !                   produced by DstrMngr.
 !      * NOTES: The variables DirichIS, CauchyIS and NeummanIS have redundant information that is 
 !             contained in the Tplgy but anyway is needed to transfer the information.

    TYPE Geometry
        DM                              :: DstrMngr
        Vec                             :: Tplgy
        IS                              :: DirichIS,CauchyIS,NeummanIS
    END TYPE Geometry


 !  - Tensor: It's a data structure that represent a full tensor, nine components, of any property.
 !      * NOTES: To use this structure is first needed fill the xx, yy, zz, xy, xz, and yz 
 !               components and then use TargetFullTensor subroutine to have (xy=yx, xz=zx, and 
 !               yz=zy). This process is required over every tensor created; thereafter, is 
 !               guaranteed the symmetry on the tensor.

    TYPE Tensor
        PetscReal                       :: xx,yy,zz,xy,xz,yz
        PetscReal,POINTER               :: yx,zx,zy
    END TYPE Tensor

 !  - Position: It's a data structure that describe a spatial position.
 !      + i,j,k: Integer that describe a global position on x, y, and z axes respectively.

    TYPE Position
        PetscInt                        :: i,j,k
    END TYPE Position

 !  - Property: It's a data structure wich describe completly a cell on the geometry.
 !      + Pstn: It's a Position structure that describe the global position on the grid.
 !      + StencilType: It's an integer used to describe the type of the stencil used on the model.
 !          0: Not defined. (Default)
 !          1: Star stencil. Stencil based on 6 neighbors.
 !          2: Partial box stencil. Stencil based on 18 neighbors. 
 !          3: Box stencil. Stencil based on 26 neighbors.
 !      + StencilTplgy: It's an Array of integers that contains an identificator of the type of the
 !                      topology on each cell of the stancil. The array is ordered numbering the 
 !                      cells of the stencil from upper to lowest layer, then from the lowest to the
 !                      highest values on the x and on the y axis respectively.
 !          0: Inactive cell.
 !          1: Active cell.
 !          2: Dirichlet boundary contition cell.
 !          3: Neumman on x boundary condition cell (dh/dx=0).
 !          4: Neumman on y boundary condition cell (dh/dy=0).
 !          5: Neumman on z boundary condition cell (dh/dz=0).
 !      + dx,dy,dz: A real that describe the size of the central cell in the directions on x, y,
 !                  and z axis.
 !      + dxB,dyB,dzB: A real that describe the size of the backward cell in the directions on x, 
 !                  y, and z axis.
 !      + dxF,dyF,dzF: A real that describe the size of the forward cell in the directions on x, 
 !                  y, and z axis.
 !      + CvtOnBlock: It's a boolean that shows if the block has the conductivity represented 
 !                    on the block.
 !      + CvtBlock: It's a Tensor of conductivities of the medium of the block.
 !      + CvtOnInterface: It's a boolean that shows if the block has the conductivity represented 
 !                        on the interfaces.
 !      + CvtBx,CvtBy,CvtBz: It's a tensor of conductivites on the interface of the block in
 !                           backward direction of the cartesian axis respectively.
 !      + CvtFx,CvtFy,CvtFz: It's a tensor of conductivites on the interface of the block in forward 
 !                           direction of the cartesian axes respectively.
 !      * NOTES: Do not use this structure to define each cell on a field of properties, it is 
 !               because this structure has a redundant data that another one already has too. 
 !               Instead, use this one as temporal structure to have everything as you need of the 
 !               cell in hand.

    TYPE Property
        TYPE(Position)                  :: Pstn
        PetscInt                        :: StnclType=0
        PetscInt,ALLOCATABLE            :: StnclTplgy(:)
        PetscReal                       :: dx,dy,dz,dxB,dxF,dyB,dyF,dzB,dzF
        PetscBool                       :: CvtOnInterface=.FALSE.,CvtOnBlock=.FALSE.
        TYPE(Tensor)                    :: CvtBlock,CvtBx,CvtFx,CvtBy,CvtFy,CvtBz,CvtFz
    END TYPE Property

 !  - ConductivityField: It's a data structure that store the field of conductivities of every block.
 !                       This structure may store conductiviies defined by field zones or where each
 !                       conductivity block has a diferent vaule.
 !      + DefinedByZones: It's a boolean describing if the ConductivityField is stored by zones or 
 !                        not. If the boolean is .TRUE., it's necessary save the information in 
 !                        CvtZones and CvtYpe, otherwise it's necessary save the information in 
 !                        xxVec,yyVec,zzVec
 !      + DefindeByInteface: It's a boolean describing if the ConductivityField stores the
 !                           information on the interface or on each block. If the boolean is
 !                           .TRUE., the CvtType or xxVec, yyVec, and zzVec need an aditional
 !                           component on each direction in the first processor.
 !      + CvtZones: It's an Array of Tenors that contain many Tensors of conductivity as defined 
 !                  zones. This variable needs to be stored in all processors.
 !      + CvtType: It's a PETSc vector that contain a identificator by each cell, the indentificator
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
 !      + Cvt: It's a ConductivityField data structure wich contain a field of conductivity.

    TYPE PropertyField
        TYPE(ConductivityField)         :: Cvt
    END TYPE PropertyField

 !  - StepBC: It's a data structure that contain a Boundary Condition set.
 !      + Diricl: It's a PETSc vector with the values of the Dirichlet cells. It must has the same
 !                size as Dirichlet indentificators defined in Gemetry data strucure.
 !      + Cauchy: It's a PETSc vector with the values of the Cauchy cells. It must has the same size
 !                as Cauchy indentificators defined in Gemetry data strucure.
 !      * NOTES: The Index Sets defined in Geometry for each boundary condition point the value of  
 !               Dirich and Cauchy vectors to the real position on every global vector crated with
 !               DstrMngr defined in Geometry.
    TYPE StepBC
        Vec                             :: Dirich,Cauchy
    END TYPE StepBC

 !  - BoundaryConditions: It's a data structure wich contains a serie Boundary Conditions to be used 
 !                        in several times or zones of time.
 !      + Step: It's an StepBc array with several values of Boundary Conditions.

    TYPE BoundaryConditions
        TYPE(StepBC),ALLOCATABLE        :: Step(:)
    END TYPE BoundaryConditions

 !  - StensilVar: It's a auxilar data structure to build the matrix.
 !      + idx_rws: It's a PETSc MatStencil structure that conitans the information of the rows that
 !                 will be added on the matrix. The allocation to the first component is 4 ever, to 
 !                 the second one dependes on how many rows are necessary to modify, usualy one.
 !      + idx_clmn: It's a PETSc MatStencil structure that conitans the information of the columns 
 !                  that will be modified on the matrix. The allocation to the first component is 4 
 !                  ever, to the second one dependes on how many colums are necessary to modify.
 !      + Values: It's an array that store the values to be added on the position indicated by 
 !                idx_rws and idx_clmns.
 !      + Size: It's an integer taht says how many values will be modified on the matrix.
 !      * NOTES: An example of usage of the StencilVar of 3 values in a row is the following:
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
 !!          StencilVar%Values(2)=2.0                    !    on the positions defined with idx_rws
 !!          StencilVar%Values(3)=1.0                    !    and idx_clmn.
 !                                  
 !               Finally, the StencilVar needs to be added to the matrix with MatSetValuesStencil
 !               that is a PETSc function.

    TYPE StencilVar
        MatStencil,ALLOCATABLE          :: idx_rws(:,:),idx_clmns(:,:)
        PetscReal,ALLOCATABLE           :: Values(:)
        PetscInt                        :: Size
    END TYPE StencilVar

 !  - RunOptionsVar: It's a data structure that contains all options related with the running.
 !      + Time: It's a boolean that define if the running is steady or transient. The default value
 !              is steady.
 !      + Scheme: It's an integer that define the scheme stencil to solve the underground flow 
 !                equation.
 !          1: Traditional solution. It's used by MODFLOW. It uses an interface conductivity defined
 !             by the armonic. Default use if is not indicated.
 !          2: Li model. It uses an interface conductivity to build the stencil.
 !          3: ANISOFLOWPACK model. -----

    TYPE RunOptionsVar
        PetscBool                       :: Time=.FALSE.
        PetscInt                        :: Scheme=1
    END TYPE RunOptionsVar

 !  - InputTypeVar: It's a collection of integer that define the type of input to be used in the
 !                 program.
 !      + Gmtry: It's a integer that define a type of file that will provide the domain size and the
 !               cell center coordenates.
 !          1: Defined by Blessent. An example is provided in ../ex/Blessent/in/tsim_USMH.asc
 !          2: Defined by Perez. An exaple is provided in ../ex/Perez/in/sanpck.domnRST
 !      + Tplgy: It's a integer that define a type of file that will provide the identificators to
 !               each cell.
 !          1: Defaul topology, it doesn't need a file. Every face of the domain is a Neumman
 !             condition but on the fisrt layer where the bpundary is dirichlet.
 !      + Cvt: It's a integer that define a type of file that will provide conductivity to each cell
 !             or zone.
 !          1: It's a file that provide block conductivities by zones. An example is provided in ../
 !             ex/Blessent/in/matrix.mprops
 !      + BC: It's

    TYPE InputTypeVar
        PetscInt                        :: Gmtry,Tplgy,Cvt,BC
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