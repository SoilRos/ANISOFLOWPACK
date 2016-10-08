MODULE ANISOFLOW_Types

! ANISOFLOW_Types it's a module that contains the basic data structures to ANISOFLOWPACK.

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscis.h>
#include <petsc/finclude/petscdm.h>


 !  - Geometry: It's a data structure that manages every information needed related to geometry.
 !    > VARIABLES: Scale, DataMngr, Tplgy, pTplgy, DirichIS, CauchyIS, NeummanIS.
 !      + Scale: It's a PETSc Integer that indentify the geometry scale
 !          1: Default scale when no other scales are being used, fine scale otherwise. Collective
 !          2: Upscale geometry. Collective
 !          3: Fine Block, usually an extraction of Fine-scale. Not Collective
 !      + DataMngr: It's a PETSc manager data for structured grid in 2 or 3 dimensions, providing general 
 !                  information on the grid to be managed in parallel programming.
 !      + Tplgy: It's a PETSc vector type produced by DataMngr that contains a topology identifier
 !               in each cell. It is changed over the time steps by the boundary condition file.
 !          0: Inactive cell.
 !          1: Active cell.
 !          2: Dirichlet boundary condition cell.
 !          3: Source in cell Q (It's treated as active cell).
 !          4: Cauchy boundary condition cell.
 !      + pTplgy: It's a PETSc vector which has the a permanent form of Tplgy; Provided by the topology file.
 !      + SizeTplgy: It's an array of integers that quantify the indentifiers of each Tplgy type decribed
 !                   where inactive cells quantifier are not counted.
 !          1: Active cell quantifier.
 !          2: Dirichlet boundary condition cell quantifier.
 !          3: Source in cell quantifier.
 !          4: Cauchy boundary condition quantifier.
 !      + x,y,z: It's a PETSc vector that contains the grid coordinates on each direction.
 !    > NOTES: The variable SizeTplgy doesn't have inactive cell identificator because it can be obtained
 !             from the others.

    TYPE Geometry
        PetscInt                        :: Scale=1
        DM                              :: DataMngr
        Vec                             :: Tplgy,pTplgy
        PetscInt                        :: SizeTplgy(4)
        Vec                             :: x,y,z
    END TYPE Geometry


 !  - Tensor: It's a data structure that represents a full tensor, nine components, of any property.
 !    > VARIABLES: xx, yy, zz, xy, xz, yz, yx, zx, zy.
 !      + xx,yy,zz,xy,xz,yz,yx,zx,zy: Property components.
 !    > NOTES: To use this structure is first needed filling the xx, yy, zz, xy, xz, and yz 
 !               components and then use TargetFullTensor subroutine to have symmetry (because is a 
 !               continuum media: xy=yx, xz=zx, and yz=zy). This process is required over every tensor 
 !               created; after that is guaranteed the symmetry on the tensor.

    TYPE Tensor
        PetscReal                       :: xx,yy,zz,xy,xz,yz
        PetscReal,POINTER               :: yx,zx,zy
    END TYPE Tensor

 !  - Position: It's a data structure that describes a spatial position.
 !      + i,j,k: It's an integer that describes a global position on x, y, and z axes respectively.

    TYPE Position
 !    > VARIABLES: i,j,k
        PetscInt                        :: i,j,k
    END TYPE Position

 !  - Property: It's a data structure that describes completely a cell and its boundaries.
 !    > VARIABLES: Pstn, StencilType, StencilTplgy, dx, dy, dz, Cvt.
 !      + Pstn: It's a Position structure that describes the global position of the cell on the grid.
 !      + StencilType: It's an integer used to describe the type of the stencil used on the model.
 !          0: Not defined. (Default)
 !          1: Star stencil. Stencil based on 6 neighbors.
 !          2: Partial box stencil. Stencil based on 18 neighbors. 
 !      + StencilTplgy: It's an array of integers that contains topology identifiers on each cell 
 !                      of the stencil. The array is ordered numbering the cells of the stencil from 
 !                      the lowest to the highest values on the x and the y axis respectively, then 
 !                      from upper to lowest layer. Note that it can be an of 6 or 18 itegers depending
 !                      on the StencilType.
 !          0: Inactive cell.
 !          1: Active cell.
 !          2: Dirichlet boundary condition cell.
 !          3: Source in cell Q (It's treated as active cell).
 !          4: Cauchy boundary condition
 !      + dx,dy,dz: It's an array of reals that describes the size of the cells directions on the x, y,
 !                  and z. It has the same order of StencilTplgy
 !      + Cvt: It's an array of Tensors of conductivities that represents the conductivities of the
 !             center of the each block. It has the same order of StencilTplgy
 !    > NOTES: Don't use this structure to define each cell on a field of properties, it is 
 !             because this structure has a redundant data that another one already has too. 
 !             Instead, use this one as a temporal structure to have everything as you need of the 
 !             cell in hand.

    TYPE Property
        TYPE(Position)                  :: Pstn
        PetscInt                        :: StnclType=0
        PetscInt,ALLOCATABLE            :: StnclTplgy(:)
        PetscReal,ALLOCATABLE           :: dx(:),dy(:),dz(:)
        TYPE(Tensor),ALLOCATABLE        :: Cvt(:)
    END TYPE Property

 !  - ConductivityField: It's a data structure that stores the field of conductivities of every
 !                       block. This structure may stores conductivities defined by field zones or
 !                       where each conductivity block has a different value.
 !    > VARIABLES: DefinedBy, ZoneID, Zone, Cell.
 !      + DefinedBy: It's an integer describing if the ConductivityField is stored by zones or 
 !                   by cell. The information stored by zones can be saved directly in ZoneID within
 !                   ConductivityField (DefinedByCvtZones) or can be pointed to the properties 
 !                   zones identificators (DefinedByPptZones).
 !      + ZoneID: It's a PETSc vector that contains an identifier for each cell, the identifier 
 !                value correspond to a zone value. If DefinedBy=2, this ZoneID is pointed to 
 !                ZoneID in PropertiesField structure. It's saved as local vector
 !      + Zone: It's an array of Tensor that contain as Tensors of conductivities as defined 
 !              zones. This variable has to be stored in all processors.
 !      + xx,yy,zz,xy,xz,yz: It's a PETSc vector that contains the each componet conductivity of each cell. 
 !                           Used when DefinedBy==3. It's saved as a local vector to have easy acces when 
 !                           the program is gathering local information (with its neighbor).
 !      + yx,zx,zy: The mirror of xy,xz,yz

    TYPE ConductivityField
        PetscInt                                :: DefinedBy=0 ! 1: DefinedByCvtZones, 2:DefinedByPptZones,3:DefinedByCell
        ! Conductivity defined by zones (Local):
        Vec                                     :: ZoneID
        TYPE(Tensor),ALLOCATABLE                :: Zone(:)
        ! Conductivity defined on every cell (Local):
        Vec                                     :: xx,yy,zz,xy,xz,yz
        Vec,POINTER                             :: yx,zx,zy
    END TYPE ConductivityField

 !  - SpecificStorageField: It's a data structure that stores the field of specific storage of every
 !                          block. This structure may stores specific storage defined by field zones or
 !                          where each specific storage block has a different value.
 !    > VARIABLES: DefinedBy, ZoneID, Zone, Cell.
 !      + DefinedBy: It's an integer describing if the SpecificStorageField is stored by zones or 
 !                   by cell. The information stored by zones can be saved directly in ZoneID within
 !                   SpecificStorageField (DefinedByStoZones) or can be pointed to the properties 
 !                   zones identificators (DefinedByPptZones).
 !      + ZoneID: It's a PETSc vector that contains an identifier for each cell, the identifier 
 !                value correspond to a zone value. If DefinedBy=2, this ZoneID is pointed to 
 !                ZoneID in PropertiesField structure. It's saved as local vector
 !      + Zone: It's a PETSc vector that contains the specific storages as defined 
 !              zones. This variable has to be stored in all processors.
 !      + Cell: It's a PETSc vector that contains the specific storage of each cell. Used when DefinedBy==3.
 !   > NOTES: The difference between Global and Local it's because Sto%ZoneID can be a pointer of Ppt%ZoneID
 !             which has to be Local, but Cell is needed in its Global form to avoid change from Local to 
 !           Global on each time step. 

    TYPE SpecificStorageField
        PetscInt                                :: DefinedBy=0 ! 1: DefinedByStoZones, 2:DefinedByPptZones,3:DefinedByCell
        ! Specific Storage defined by zones (Local):
        Vec                                     :: ZoneID
        Vec                                     :: Zone
        ! Specific Storage defined on every cell (Global).:
        Vec                                     :: Cell
    END TYPE SpecificStorageField

 !  - ProperyField: It's a data structure wich contains fields of differents properties.
 !    > VARIABLES: Cvt,Sto, DefinedByPptZones
 !      + Cvt: It's a ConductivityField data structure wich contain a field of conductivities.
 !      + Sto: It's a StorageField data structure wich contain a field of specific storages.
 !      + DefinedByPptZones: It's a Boolean which says if the property field is stored by zones;
 !                           in such case, the variable ZoneID will characterize each cell by zones
 !                           This variable helps to the user to avoid defining a field zone to every
 !                           single property.
 !      + ZoneID: It's a PETSc vector that contains an identifier for each cell, the identifier 
 !                value correspond to a zone value.
 !    > NOTES: In some cases the Cvt%ZoneID and/or Sto%ZoneID will point to ZoneID of the PropertyField.
 !             Depends on how the user wants to input his variables

    TYPE PropertiesField
        TYPE(ConductivityField)         :: Cvt
        TYPE(SpecificStorageField)      :: Sto
        ! Property defined by zones (Local):
        PetscBool                       :: DefinedByPptZones=.FALSE.
        Vec                             :: ZoneID
    END TYPE PropertiesField

 !  - TimeZoneVar: It's a data structure a zone of time to be modeled.
 !    > VARIABLES: SizeTime,Time
 !      + SizeTime: It's an Integer that says the amount of time values.
 !      + Time: It's an Array that contains the time values to be used in this Time Zone.
 !    > NOTES: The TimeZoneVar represents one zone, such zone can has several time steps.
 !             It's characterized because the boundary condition never changes within the zone.

    TYPE TimeZoneVar
        PetscInt                        :: SizeTime
        PetscReal,ALLOCATABLE           :: Time(:)
    END TYPE TimeZoneVar

 !  - BoundaryConditions: It's a data structure wich contains a series Boundary Conditions to be  
 !                        used in several zones of time.
 !    > VARIABLES: SizeTimeZone, SizeDirich, SizeSource, SizeCauchy, Dirich, Source, Cauchy,
 !                 TimeZone.
 !      + SizeTimeZone: It's an Integer that says the amount of time zones.
 !      + SizeDirich: It's an array of integers that says the size of Dirich PETSc vector for each time zone.
 !      + SizeSource: It's an array of integers that says the size of Neumman PETSc vector for each time zone.
 !      + SizeCauchy: It's an array of integers that says the size of Cauchy PETSc vector for each time zone.
 !      + Dirich: It's a PETSc vector with the values of the Dirichlet cells to each time zone.
 !      + Neumman: It's a PETSc vector with the values of the Neumman cells to each time zone.
 !      + CauchyC: It's a PETSc vector with the values of C part of the Cauchy cells to each time zone.
 !      + CauchyHe: It's a PETSc vector with the values of He part the Cauchy cells to each time zone.
 !                  (Saved in its negative form to apply the vector directly to the system.)
 !      + DirichIS: It's a PETSc index set that has a map between Dirichlet information and vecs 
 !                  produced by DataMngr.
 !      + SourceIS: It's a PETSc index set that has a map between Source information and vecs 
 !                   roduced by DataMngr.
 !      + CauchyIS: It's a PETSc index set that has a map between Cauchy information and vecs 
 !                  produced by DataMngr.
 !    > NOTES: The Index Sets defined in Geometry for each boundary condition. And Points the value
 !             of Dirich, Neumman, and Cauchy vectors to the real position of every global vector 
 !             created with DataMngr defined in Geometry.

    TYPE BoundaryConditions
        PetscInt                        :: SizeTimeZone
        PetscInt,ALLOCATABLE            :: SizeDirich(:),SizeSource(:),SizeCauchy(:)
        Vec,ALLOCATABLE                 :: Dirich(:),Source(:),CauchyC(:),CauchyHe(:)
        IS,ALLOCATABLE                  :: DirichIS(:),SourceIS(:),CauchyIS(:)
        TYPE(TimeZoneVar),ALLOCATABLE   :: TimeZone(:)
    END TYPE BoundaryConditions

 !  - StensilVar: It's a auxilar data structure to build the matrix.
 !    > VARIABLES: idx_rws, idx_clmn, idx_val, Size.
 !      + idx_rws: It's a PETSc MatStencil structure that contains the information of the rows that 
 !                 will be modified in the matrix. The allocation to the first component is always 
 !                 four, to the second one it depends on how many rows are necessary to change, 
 !                 usually one.
 !      + idx_clmn: It's a PETSc MatStencil structure that contanis the information of the columns 
 !                  that will be modified on the matrix. The allocation to the first component is  
 !                  always four, to the second one it depends on how many colums are necessary to 
 !                  change.
 !      + idx_val: It's an array that stores the values to be added to the position indicated by 
 !                 idx_rws and idx_clmns.
 !      + idx_size: It's an integer that says how many values will be modified on the matrix.
 !    > NOTES: An example of usage of the StencilVar of 3 values in a row is the following:
 !          
 !!          StencilVar%idx_size=3                           ! Amount of values of the Stencil
 !!
 !!          ALLOCATE(StencilVar%idx_rws(4,1))           ! One row to modify.
 !!          ALLOCATE(StencilVar%idx_clmn(4,3))          ! Four values of the row to modify.
 !!          ALLOCATE(StencilVar%idx_val(StencilVar%idx_size))! Amount of values of the Stencil
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
 !!          StencilVar%idx_val(1)=1.D0                   ! Values of the stencil that will be added 
 !!          StencilVar%idx_val(2)=2.D0                   !    to the positions defined with idx_rws
 !!          StencilVar%idx_val(3)=1.D0                   !    and idx_clmn.
 !                                  
 !             Finally, the StencilVar needs to be added to the matrix with MatSetValuesStencil
 !             that is a PETSc function.

    TYPE StencilVar
        MatStencil,ALLOCATABLE          :: idx_rws(:,:),idx_clmns(:,:)
        PetscReal,ALLOCATABLE           :: idx_val(:)
        PetscInt                        :: idx_size
    END TYPE StencilVar

 !  - RunOptionsVar: It's a data structure that contains all options related to the program execution.
 !    > VARIABLES: Time, Scheme.
 !      + Time: It's a boolean that defines if the excecution is steady or transient. The default value
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
 !                  program.
 !    > VARIABLES: Gmtry, Tplgy, Cvt, BC.
 !      + Gmtry: It's an integer that defines a type of file that will provide the domain size and
 !               the cell center coordinates.
 !          1: Defined by Blessent. An example is provided in "***"
 !          2: Defined by Perez. An example is provided in "../ex/Perez/in/sanpck.domnRST"
 !      + Tplgy: It's an integer that defines a type of file that will provide the identifiers to
 !               each cell.
 !          1: Default. From a file which has a list of integers with the Tplgy identifiers ordered 
 !             the lowest to the highest values on the x and the y axis respectively, then from
 !             upper to lowest layer.
 !          2: This type doesn't need an input file; the border of first layer is dirichlet otherwise is active.
 !          3: This type doesn't need an input file; All borders are dirichlet otherwise is active.
 !      + Cvt: It's an integer that defines a type of file that will provide conductivity to each 
 !             cell or zone.
 !          1: The Cvt is defided by a pair of files that provide block conductivities by zones and 
 !             that provide an identifier of zone to each cell. An example is provided in "../ex/
 !             Blessent/in/matrix.mprops" and "***". The first file can 
 !             be entered to describe the conductivities as well as all properties.
 !          2: The Cvt is defined by one file that contains a list of conductivities of each block.
 !             It has the same order as of Tplgy input.
 !      + BC: It's a integer that define a type of file that will provide the boundary conditions.
 !          1: It's a file that only provide boundary conditions and their position on the grid. An
 !             example is provided in "***". In this input type the potition have to have a consecutive 
 !             ordering, according to Tplgy input.
 !          2: It's a file that only provide dirichlet condition and their position on the grid. An
 !             example is provided in "***". The potition have to be described by its i, j and k position
 !             on the grid.

    TYPE InputTypeVar
        PetscInt                        :: Gmtry,Tplgy,Cvt,Sto,BC,InitSol
    END TYPE InputTypeVar

 !  - OutputTypeVar: It's a collection of integer that defines the type of ouput to be used in the
 !                  program.
 !    > VARIABLES: Sol.
 !      + Sol: It's an integer that defines a type of ouput to be used in solution.
 !          1: Binary
 !          2: ASCII
 !          3: HDF5

    TYPE OutputTypeVar
        PetscInt                        :: Sol=3,Tplgy=0,Cvt=0,Sto=0,Ppt=0
    END TYPE OutputTypeVar

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

    END FUNCTION TargMirrorValue

!!   - TargMirrorValue: It's a auxilar function to target a real that will be pointed later.
!!     > IN: Value
!!       + Values: It's a real.
!!     > OUT: Value
!!       + Values: It's the same input real but now as a target.

    FUNCTION TargPetscVec(x)

        IMPLICIT NONE

        Vec,INTENT(IN),TARGET     :: x
        Vec,POINTER               :: TargPetscVec
        TargPetscVec => x

    END FUNCTION TargPetscVec

END MODULE ANISOFLOW_Types









