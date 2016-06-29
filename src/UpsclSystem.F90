MODULE ANISOFLOW_UpsclSystem

    IMPLICIT NONE

CONTAINS

!This subrutine manage all subrutines related with upscaling procedures

SUBROUTINE GetUpsclSystem(FineGmtry,PptFld,BCFld,ierr)
    USE ANISOFLOW_Types
    USE ANISOFLOW_Geometry
    USE ANISOFLOW_BuildSystem
    USE ANISOFLOW_Solver

#include <petsc/finclude/petscsys.h>

    IMPLICIT NONE

    PetscErrorCode,INTENT(INOUT)               :: ierr
    TYPE(Geometry),INTENT(INOUT)               :: FineGmtry
    TYPE(PropertyField),INTENT(INOUT)          :: PptFld
    TYPE(BoundaryConditions),INTENT(INOUT)     :: BCFld

    TYPE(RunOptionsVar)                        :: RunOptions
    TYPE(Geometry)                             :: UpsclGmtry,FineBlockGmtry
    PetscInt                                   :: u,i,j,k,w,Upsclwidth(3),skin




    !TODO: Get geometry at coarse scale, arreglar la interface para esto
    CALL GetUpsclGeometry(UpsclGmtry,ierr) !TODO: Revisar que las mallas x,y,z cuadren entre si
    ! Operaciones con el UpsclGmtry

    CALL DMDAGetCorners(UpsclGmtry%DataMngr,UpsclCorn(1),UpsclCorn(2),UpsclCorn(3),UpsclWidthL(1),      &
            & UpsclWidthL(2),UpsclWidthL(3),ierr)

    ! recorre la malla gruesa en cada procesador
    DO k=UpsclCorn(3),UpsclCorn(3)+UpsclWidthL(3)-1
        DO j=UpsclCorn(2),UpsclCorn(2)+UpsclWidthL(2)-1
            DO i=UpsclCorn(1),UpsclCorn(1)+UpsclWidthL(1)-1

                CALL GetFineBlockGmtry(FineGmtry,UpsclGmtry,i,j,k,FineBlockGmtry,ierr)
!                 DO w=1,8
!                     CALL GetSystem(FineGmtry,PptFld,BCFld,i,j,A,b,x,ierr)
!                     CALL SolveSystem(FineGmtry,PptFld,BCFld,A,b,x,ierr)
!                 END DO
            END DO
        END DO
    END DO

!     FineGmtry=UpsclGmtry



END SUBROUTINE GetUpsclSystem

!This subrutine build the fine geometry related to a block of coarse geometry
!i=Bloque iesimo
!j=Bloque siguiente en direccion Y
!Bloque siguiente en direccion Z

SUBROUTINE GetFineBlockGmtry(FineGmtry,UpsclGmtry,Ui,Uj,Uk,FineBlockGmtry,ierr)

    USE ANISOFLOW_Types

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
    
    IMPLICIT NONE

    PetscErrorCode,INTENT(INOUT)                :: ierr
    TYPE(Geometry),INTENT(IN)                   :: FineGmtry,UpsclGmtry
    TYPE(Geometry),INTENT(OUT)                  :: FineBlockGmtry
    PetscInt,INTENT(IN)                         :: Ui,Uj,Uk !Upscale indices

    PetscInt                                    :: FiB,FiF,FjB,FjF,FkB,FkF !Fine indices backward and forward
    PetscInt                                    :: Bi,Bj,Bk !Block indices
    PetscInt                                    :: GridSize(3),Width(3)
    PetscReal                                   :: xB,xF,yB,yF,zB,zF
    Vec                                         :: TmpTplgy
    PetscReal,POINTER                           :: TmpTplgyArray(:,:,:)

    ! Skin used in upscaling procedure
    PARAMETER(Skin=5) !TODO: Interface GetUpsclSkin

    ! Get grid -----------------------------
!     xB=UpsclGmtry%x(Ui)
!     xF=UpsclGmtry%x(Ui+1)

!     yB=UpsclGmtry%y(Uj)
!     yF=UpsclGmtry%y(Uj+1)

!     zB=UpsclGmtry%z(Uk)
!     zF=UpsclGmtry%z(Uk+1)


    FiB=BinarySearch(FineGmtry%x,xB)
    FiF=BinarySearch(FineGmtry%x,xF)

    GridSize(1)=FiF-FiB+2*Skin
    ALLOCATE(FineBlockGmtry%x(GridSize(1)))
    FineBlockGmtry%x=FineGmtry%x(FiB-Skin:FiF+Skin)


    FjB=BinarySearch(FineGmtry%y,yB)
    FjF=BinarySearch(FineGmtry%y,yF)

    GridSize(2)=FjF-FjB+2*Skin
    ALLOCATE(FineBlockGmtry%y(GridSize(2)))
    FineBlockGmtry%y=FineGmtry%y(FjB-Skin:FjF+Skin)

    ! TODO: Condicional para 2D
    FkB=BinarySearch(FineGmtry%z,zB)
    FkF=BinarySearch(FineGmtry%z,zF)

    GridSize(3)=FkF-FkB+2*Skin
    ALLOCATE(FineBlockGmtry%z(GridSize(3)))
    FineBlockGmtry%z=FineGmtry%z(FkB-Skin:FkF+Skin)


    ! Get DMDA ------------------------------
    Width=GridSize-1
    ! TODO: Condicional 2D
    CALL DMDACreate3d(PETSC_COMM_SELF,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,      &
        & DM_BOUNDARY_GHOSTED,Stencil,Width(1),Width(2),Width(3),     &
        & PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,PETSC_NULL_INTEGER,       &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,FineBlockGmtry%DataMngr,ierr)

    ! Get Topology -------------------------

    ! Get topology with every border as dirichlet
    CALL GetTopology_3(PETSC_COMM_SELF,FineBlockGmtry%DataMngr,FineBlockGmtry%Tplgy,FineBlockGmtry%DirichIS,FineBlockGmtry%SourceIS,FineBlockGmtry%CauchyIS,ierr)


END SUBROUTINE GetFineBlockGmtry

! SUBROUTINE GetUpscldBC(Gmtry,PptFld,BCFld,ierr)


! END SUBROUTINE GetUpscldBC


! TODO: Test recursive function vs while function

RECUSRIVE FUNCTION BinarySearch(Array,Goal) RESULT(Index)

#include <petsc/finclude/petscsys.h>

    IMPLICIT NONE

    PetscReal,DIMENSION(:),INTENT(IN)   :: Array ! Ordered array!!
    PetscReal,INTENT(IN)                :: Goal
    PetscInt,INTENT(OUT)                :: Index

    PetscInt                            :: Middle

    IF (SIZE(Array).EQ.1).AND.(Array(1).NE.Goal) THEN
        PRINT, "Error in BinarySearch: Goal don't exist"
        STOP
    END IF

    Middle=(1+SIZE(Array))/2
    IF (Goal.EQ.Array(Middle)) THEN
        Index=Middle
        RETURN
    ELSE IF (Goal.GT.Array(Middle)) THEN
        Index=BinarySearch(Array(1:Middle-1),Goal)
        RETURN
    ELSE
        Index=BinarySearch(Array(Middle+1:),Goal)
        Index=Index+Middle
        RETURN
    END IF

END FUNCTION BinarySearch

END MODULE ANISOFLOW_SetUpsclSystem

