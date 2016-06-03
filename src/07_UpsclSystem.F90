MODULE ANISOFLOW_UpsclSystem

    IMPLICIT NONE

CONTAINS

!This subrutine manage all subrutines related with upscaling procedures
!

SUBROUTINE GetUpsclSystem(Gmtry,PptFld,BCFld,ierr)
    USE ANISOFLOW_Types
    USE ANISOFLOW_Geometry
    USE ANISOFLOW_BuildSystem
    USE ANISOFLOW_Solver

#include <petsc/finclude/petscsys.h>


    PetscErrorCode,INTENT(INOUT)               :: ierr
    TYPE(Geometry),INTENT(INOUT)               :: Gmtry
    TYPE(PropertyField),INTENT(INOUT)          :: PptFld
    TYPE(BoundaryConditions),INTENT(INOUT)     :: BCFld

    TYPE(RunOptionsVar)                        :: RunOptions
    TYPE(Geometry)                             :: UpsclGmtry,FineBlockGmtry
    PetscInt                                   :: u,i,j,k,w,Upsclwidth(3),skin

    !Skin used in upscaling procedure
    PARAMETER(skin=5)

    ! It obtains run options.
    CALL GetRunOptions(RunOptions,ierr)

    IF (RunOptions%Scheme.EQ.3) THEN
        !Get geometry at coarse scale
        CALL GetGeometry(UpsclGmtry,ierr)
        ! Operaciones con el UpsclGmtry
        CALL DMDAGetInfo(UpsclGmtry%DataMngr,PETSC_NULL_INTEGER,Upsclwidth(1),Upsclwidth(2),&
        & Upsclwidth(3),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,                 &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,        &
        & PETSC_NULL_INTEGER,ierr)

        DO i=1,Upsclwidth(1)*Upsclwidth(2)*Upsclwidth(3)
            j=i+Upsclwidth(1)
            k=i+Upsclwidth(1)*Upsclwidth(2)
            !i,j,k coordenada
            CALL GetFineBlockGmtry(Gmtry,UpsclGmtry,i,j,k,FineBlockGmtry,ierr)
            DO w=1,8
                CALL GetSystem(Gmtry,PptFld,BCFld,i,j,A,b,x,ierr)
                CALL SolveSystem(Gmtry,PptFld,BCFld,A,b,x,ierr)
            END DO
        END DO

        Gmtry=UpsclGmtry

    END IF


END SUBROUTINE GetUpsclSystem

!This subrutine construct the fine geometry related to a block of coarse geometry
!i=Bloque iesimo
!j=Bloque siguiente en direccion Y
!Bloque siguiente en direccion Z

SUBROUTINE GetFineBlockGmtry(Gmtry,UpsclGmtry,i,j,k,skin,FineBlockGmtry,ierr)

    USE ANISOFLOW_Types

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)               :: ierr
    TYPE(Geometry),INTENT(IN)                  :: Gmtry,UpsclGmtry
    TYPE(Geometry),INTENT(OUT)                 :: FineBlockGmtry
    PetscInt,INTENT(IN)                        :: i,j,k,skin

    Vec                                        :: x,y,z
    TYPE(Position)                             :: a,b,c,d
    PetscInt                                   :: sizea,aux


    ! Array dimension in x direction
    size=b-a+2*skin
    CALL VecCreateSeq(PETSC_COMM_SELF,size,x,ierr)
    ! Array dimension in y direction
    size=c-a+2*skin
    CALL VecCreateSeq(PETSC_COMM_SELF,size,y,ierr)
    ! Array dimension in z direction
    size=d-a+2*skin
    CALL VecCreateSeq(PETSC_COMM_SELF,size,z,ierr)

    !a centro, b-x, c-y, d-z
    aux=1
    DO i=a-skin,b+skin
        x(aux)=Gmtry%x(i)
        aux=aux+1
    END DO

    aux=1
    DO i=a-skin,c+skin
        y(aux)=Gmtry%y(i)
        aux=aux+1
    END DO

    aux=1
    DO i=a-skin,d+skin
        z(aux)=Gmtry%x(i)
        aux=aux+1
    END DO

    UpsclGmtry%x(i)
    UpsclGmtry%x(j)
    UpsclGmtry%x(k)

    FineBlockGmtry%DataMngr=
    FineBlockGmtry%Tplgy=
    FineBlockGmtry%x=x
    FineBlockGmtry%y=y
    FineBlockGmtry%z=z
    FineBlockGmtry%DirichIS=
    FineBlockGmtry%SourceIS=
    FineBlockGmtry%CauchyIS=

END SUBROUTINE GetFineBlockGmtry

SUBROUTINE GetUpscldBC(Gmtry,PptFld,BCFld,ierr)


END SUBROUTINE GetUpscldBC

END MODULE ANISOFLOW_SetUpsclSystem

