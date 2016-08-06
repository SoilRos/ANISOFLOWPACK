MODULE ANISOFLOW_BuildSystem

    IMPLICIT NONE

CONTAINS

SUBROUTINE GetSystem(Gmtry,PptFld,BCFld,TimeZone,TimeStep,A,b,x,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,PropertiesField,BoundaryConditions,RunOptionsVar
    USE ANISOFLOW_Geometry, ONLY : UpdateGmtry
    USE ANISOFLOW_Interface, ONLY : GetRunOptions,GetVerbose

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Geometry),INTENT(INOUT)            :: Gmtry
    TYPE(PropertiesField),INTENT(IN)        :: PptFld
    TYPE(BoundaryConditions),INTENT(IN)     :: BCFld
    PetscInt,INTENT(IN)                     :: TimeZone,TimeStep
    Mat,INTENT(INOUT)                       :: A
    Vec,INTENT(INOUT)                       :: b,x

    TYPE(RunOptionsVar)                     :: RunOptions
    CHARACTER(LEN=200)                      :: EventName,ClassName
    PetscBool                               :: Verbose
    PetscLogEvent                           :: Event
    PetscClassId                            :: ClassID
    PetscLogDouble                          :: EventFlops=0.d0

    ClassName="System"
    CALL PetscClassIdRegister(ClassName,ClassID,ierr)
    EventName="GetSystem"
    CALL PetscLogEventRegister(EventName,ClassID,Event,ierr)
    CALL PetscLogEventBegin(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Inizialited\n",ierr)
    
    CALL GetRunOptions(RunOptions,ierr)

    CALL BuildSystem(Gmtry,PptFld,A,ierr)
!     IF ((TimeZone.EQ.1).AND.(TimeStep.EQ.1)) THEN
!         CALL UpdateGmtry(Gmtry,BCFld%DirichIS(TimeZone),BCFld%SourceIS(TimeZone),BCFld%CauchyIS(TimeZone),ierr)
!         CALL BuildSystem(Gmtry,PptFld,A,ierr)
!         CALL GetInitSol(Gmtry,x,ierr)
! !         CALL DMCreateGlobalVector(Gmtry%DataMngr,x,ierr) ! inicializar x
!         CALL DMCreateGlobalVector(Gmtry%DataMngr,b,ierr)
!     END IF

!     IF ((RunOptions%Time).AND..NOT.((TimeZone.EQ.1).AND.(TimeStep.EQ.1))) THEN
! !         CALL UpdateGmtry(Gmtry,DirichIS,SourceIS,CauchyIS,ierr)
! !         CALL UpdateSystem(Gmtry,PptFld,A,ierr) UpdateSystem
!         CALL ApplyTimeDiff(PptFld,BCFld,TimeZone,TimeStep,A,b,x,ierr)
!         CALL ApplyDirichlet(BCFld,TimeZone,b,ierr)
!         CALL ApplySource(BCFld,TimeZone,b,ierr)
!         CALL ApplyCauchy(BCFld,TimeZone,A,b,ierr)
!     END IF


    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Finalized\n",ierr)
    
    CALL PetscLogFlops(EventFlops,ierr)
    CALL PetscLogEventEnd(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

END SUBROUTINE GetSystem

SUBROUTINE BuildSystem(Gmtry,PptFld,A,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,PropertiesField,Property,StencilVar
    USE ANISOFLOW_Properties, ONLY : GetLocalProperty
    USE ANISOFLOW_Interface, ONLY : GetVerbose

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Geometry),INTENT(IN)               :: Gmtry
    TYPE(PropertiesField),INTENT(IN)        :: PptFld
    Mat,INTENT(OUT)                         :: A

    PetscInt                                :: i,j,k,corn(3),widthL(3)
    TYPE(Property)                          :: Ppt
    TYPE(StencilVar)                        :: Stencil
    CHARACTER(LEN=200)                      :: EventName,ClassName
    PetscBool                               :: Verbose
    PetscLogEvent                           :: Event
    PetscClassId                            :: ClassID
    PetscLogDouble                          :: EventFlops=0.d0

    ClassName="System"
    CALL PetscClassIdRegister(ClassName,ClassID,ierr)
    EventName="BuildSystem"
    CALL PetscLogEventRegister(EventName,ClassID,Event,ierr)
    CALL PetscLogEventBegin(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Inizialited\n",ierr)
    

    CALL DMCreateMatrix(Gmtry%DataMngr,A,ierr)

    CALL DMDAGetCorners(Gmtry%DataMngr,corn(1),corn(2),corn(3),widthL(1),      &
            & widthL(2),widthL(3),ierr)

    DO k=corn(3),corn(3)+widthL(3)-1
        DO j=corn(2),corn(2)+widthL(2)-1
            DO i=corn(1),corn(1)+widthL(1)-1
                CALL GetLocalProperty(Gmtry,PptFld,Ppt,i,j,k,ierr)
                CALL GetStencil(Ppt,Stencil,ierr)
                
                CALL MatSetValuesStencil(A,1,Stencil%idx_rws,Stencil%Size,     &
                    & Stencil%idx_clmns,Stencil%Values,ADD_VALUES,ierr)
            END DO
        END DO
    END DO

    CALL MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)  

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"["//ADJUSTL(TRIM(EventName))//" Event] Finalized\n",ierr)
    
    CALL PetscLogFlops(EventFlops,ierr)
    CALL PetscLogEventEnd(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

END SUBROUTINE BuildSystem

SUBROUTINE GetInitSol(Gmtry,x,ierr)

    USE ANISOFLOW_Types, ONLY : InputTypeVar,Geometry
    USE ANISOFLOW_Interface, ONLY : GetInputDir,GetInputTypeInitSol,GetInputFileInitSol

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscviewer.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Geometry),INTENT(IN)               :: Gmtry
    Vec,INTENT(OUT)                         :: x

    PetscBool                   :: InputFileInitSolFlg
    TYPE(InputTypeVar)          :: InputType
    CHARACTER(LEN=200)          :: Route,InputDir,InputFileInitSol
    PetscViewer                 :: Viewer
   
    CALL GetInputDir(InputDir,ierr)
    CALL GetInputFileInitSol(InputFileInitSol,InputFileInitSolFlg,ierr)

    CALL DMCreateGlobalVector(Gmtry%DataMngr,x,ierr)

    IF (InputFileInitSolFlg) THEN
        CALL GetInputTypeInitSol(InputType,ierr)
        Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFileInitSol))
        CALL PetscObjectSetName(x,"Solution",ierr)
        IF (InputType%InitSol.EQ.1) THEN
            CALL PetscViewerBinaryOpen(PETSC_COMM_WORLD,Route,FILE_MODE_READ,Viewer,ierr)
            CALL VecLoad(x,Viewer,ierr)
            CALL PetscViewerDestroy(Viewer,ierr)
        ELSEIF (InputType%InitSol.EQ.2) THEN
            PRINT*,"ERROR: Initial solution is not able to be opend from ASCII format, you should use binary or HDF5 formats. Initial solution was set to 0."
        ELSEIF (InputType%InitSol.EQ.3) THEN
            CALL PetscViewerHDF5Open(PETSC_COMM_WORLD,Route,FILE_MODE_READ,Viewer,ierr)
            CALL VecLoad(x,Viewer,ierr)
            CALL PetscViewerDestroy(Viewer,ierr)
            PRINT*,"Inital solution ",InputFileInitSol," has been implemented properly."
        ELSE 
            PRINT*,"ERROR: InitSol is not valid. Initial solution was set to 0."
        END IF
    ELSE
        PRINT*,"WARING: Initial solution was set to 0."
    END IF

END SUBROUTINE GetInitSol

SUBROUTINE GetStencil(Ppt,Stencil,ierr)

    USE ANISOFLOW_Types, ONLY : Property,StencilVar,RunOptionsVar
    USE ANISOFLOW_Interface, ONLY : GetRunOptions

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Property),INTENT(IN)               :: Ppt
    TYPE(StencilVar),INTENT(OUT)            :: Stencil
    
    TYPE(RunOptionsVar)                     :: RunOptions
    CHARACTER(LEN=200)                      :: EventName,ClassName
    PetscLogEvent                           :: Event
    PetscClassId                            :: ClassID
    PetscLogDouble                          :: EventFlops=0.d0

    ClassName="System"
    CALL PetscClassIdRegister(ClassName,ClassID,ierr)
    EventName="GetStencil"
    CALL PetscLogEventRegister(EventName,ClassID,Event,ierr)
    CALL PetscLogEventBegin(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

    CALL GetRunOptions(RunOptions,ierr)

    ALLOCATE(Stencil%idx_rws(4,1))

    IF (RunOptions%Scheme.EQ.1) THEN
        CALL GetTraditionalStencil(Ppt,Stencil,ierr)
    ELSEIF (RunOptions%Scheme.EQ.2) THEN
        CALL GetLiStencil(Ppt,Stencil,ierr)
    ! ELSEIF (RunOptions%Scheme.EQ.3) THEN
    !     CALL AnisoflowStencil(Ppt,Stencil,ierr)
    ELSE
    
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
            & "ERROR: Run_options_scheme command must be an integer between 1 and 2\n",ierr)
        STOP
    
    END IF

    CALL PetscLogFlops(EventFlops,ierr)
    CALL PetscLogEventEnd(Event,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)


END SUBROUTINE GetStencil

!FALTANTE: Subrutina que haga calculo de matriz de almacenamiento especifico

SUBROUTINE GetTraditionalStencil(Ppt,Stencil,ierr)

    USE ANISOFLOW_Types, ONLY : Property,StencilVar

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Property),INTENT(IN)               :: Ppt
    TYPE(StencilVar),INTENT(INOUT)          :: Stencil

    PetscInt                                :: i,j,k
    PetscReal                               :: one=1.d0,zero=0.d0

    ALLOCATE(Stencil%idx_clmns(4,7))
    ALLOCATE(Stencil%Values(7))
    Stencil%Size=7

    ! It gets the position of the cell.
    i=Ppt%Pstn%i
    j=Ppt%Pstn%j
    k=Ppt%Pstn%k

    ! Rows to modify
    Stencil%idx_rws(MatStencil_i,1) = i
    Stencil%idx_rws(MatStencil_j,1) = j
    Stencil%idx_rws(MatStencil_k,1) = k
    ! Columns to modify
    Stencil%idx_clmns(MatStencil_i,:) = i
    Stencil%idx_clmns(MatStencil_j,:) = j
    Stencil%idx_clmns(MatStencil_k,:) = k

    ! Initial stencil vaulues
    Stencil%Values(:)=zero

    ! If the current cell is an active, source, or cauchy cell:
    IF ((Ppt%StnclTplgy(4).EQ.1).OR.(Ppt%StnclTplgy(4).EQ.3).OR.(Ppt%StnclTplgy(4).EQ.4)) THEN
        Stencil%idx_clmns(MatStencil_i,1) = i
        Stencil%idx_clmns(MatStencil_j,1) = j
        Stencil%idx_clmns(MatStencil_k,1) = k-1
        Stencil%Values(1)=2*Ppt%dx*Ppt%dy*Ppt%CvtBz%zz/(Ppt%dzB+Ppt%dz)

        Stencil%idx_clmns(MatStencil_i,2) = i
        Stencil%idx_clmns(MatStencil_j,2) = j-1
        Stencil%idx_clmns(MatStencil_k,2) = k
        Stencil%Values(2)=2*Ppt%dx*Ppt%dz*Ppt%CvtBy%yy/(Ppt%dyB+Ppt%dy) 

        Stencil%idx_clmns(MatStencil_i,3) = i-1
        Stencil%idx_clmns(MatStencil_j,3) = j
        Stencil%idx_clmns(MatStencil_k,3) = k
        Stencil%Values(3)=2*Ppt%dy*Ppt%dz*Ppt%CvtBx%xx/(Ppt%dxB+Ppt%dx) 

        Stencil%idx_clmns(MatStencil_i,4) = i
        Stencil%idx_clmns(MatStencil_j,4) = j
        Stencil%idx_clmns(MatStencil_k,4) = k
        Stencil%Values(4)=-2*(Ppt%dy*Ppt%dz*Ppt%CvtFx%xx/(Ppt%dxF+Ppt%dx)+Ppt%dy*Ppt%dz*Ppt%CvtBx%xx/(Ppt%dxB+Ppt%dx)  &
            & +Ppt%dx*Ppt%dz*Ppt%CvtFy%yy/(Ppt%dyF+Ppt%dy)+Ppt%dx*Ppt%dz*Ppt%CvtBy%yy/(Ppt%dyB+Ppt%dy))
        IF ((Ppt%dzF.NE.zero).AND.(Ppt%dzB.NE.zero)) THEN ! Non-NaN checker. If is not NaN, then its a 3D problem
            Stencil%Values(4)=Stencil%Values(4)-2*(Ppt%dx*Ppt%dy*Ppt%CvtFz%zz/(Ppt%dzF+Ppt%dz)+Ppt%dx*Ppt%dy*Ppt%CvtBz%zz/(Ppt%dzB+Ppt%dz))
        END IF


        Stencil%idx_clmns(MatStencil_i,5) = i+1
        Stencil%idx_clmns(MatStencil_j,5) = j
        Stencil%idx_clmns(MatStencil_k,5) = k
        Stencil%Values(5)=2*Ppt%dy*Ppt%dz*Ppt%CvtFx%xx/(Ppt%dxF+Ppt%dx)

        Stencil%idx_clmns(MatStencil_i,6) = i
        Stencil%idx_clmns(MatStencil_j,6) = j+1
        Stencil%idx_clmns(MatStencil_k,6) = k
        Stencil%Values(6)=2*Ppt%dx*Ppt%dz*Ppt%CvtFy%yy/(Ppt%dyF+Ppt%dy)

        Stencil%idx_clmns(MatStencil_i,7) = i
        Stencil%idx_clmns(MatStencil_j,7) = j
        Stencil%idx_clmns(MatStencil_k,7) = k+1
        Stencil%Values(7)=2*Ppt%dx*Ppt%dy*Ppt%CvtFz%zz/(Ppt%dzF+Ppt%dz)

    ELSEIF ((Ppt%StnclTplgy(4).EQ.2).OR.(Ppt%StnclTplgy(4).EQ.0)) THEN ! Dirichlet, Cauchy or inactive cell
        Stencil%Values(4)=-one
    END IF


END SUBROUTINE GetTraditionalStencil

 !  - GetLiStencil: It's a routine that provides a Li stencil for a cell stored in a StensilVar.
 !    > IN: Ppt.
 !      + Ppt: It's a Property structure that contains every information needed to build the stencil
 !             on that cell.
 !    > OUT: Stencil, ierr.
 !      + Stencil: It's a StencilVar data structure that contains the information to be added to the
 !                 matrix. See StencilVar for more information.
 !      + ierr: It's an integer that indicates whether an error has occurred during the call.

SUBROUTINE GetLiStencil(Ppt,Stencil,ierr)

    USE ANISOFLOW_Types, ONLY : Property,StencilVar

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Property),INTENT(IN)               :: Ppt
    TYPE(StencilVar),INTENT(INOUT)          :: Stencil

    PetscInt                                :: i,j,k
    PetscReal                               :: one=1.d0,zero=0.d0

    ! The Li stencil is based on 19 cells.
    ALLOCATE(Stencil%idx_clmns(4,19))
    ALLOCATE(Stencil%Values(19))
    Stencil%Size=19

    ! It gets the position of the cell.
    i=Ppt%Pstn%i
    j=Ppt%Pstn%j
    k=Ppt%Pstn%k

    ! Rows to modify
    Stencil%idx_rws(MatStencil_i,1) = i
    Stencil%idx_rws(MatStencil_j,1) = j
    Stencil%idx_rws(MatStencil_k,1) = k
    ! Columns to modify
    Stencil%idx_clmns(MatStencil_i,:) = i
    Stencil%idx_clmns(MatStencil_j,:) = j
    Stencil%idx_clmns(MatStencil_k,:) = k

    ! Initial stencil vaulues
    Stencil%Values(:)=zero

    ! If the current cell is an active, source, or cauchy cell:
    IF ((Ppt%StnclTplgy(10).EQ.1).OR.(Ppt%StnclTplgy(10).EQ.3).OR.(Ppt%StnclTplgy(10).EQ.4)) THEN

    ! 1-S Bloque centro-detras-superior
        ! It sets the column position on the matrix
        Stencil%idx_clmns(MatStencil_i,1) = i
        Stencil%idx_clmns(MatStencil_j,1) = j-1
        Stencil%idx_clmns(MatStencil_k,1) = k-1
        Stencil%Values(1)=Ppt%dy*Ppt%dz*Ppt%CvtBy%yz/(Ppt%dzF+2*Ppt%dz+Ppt%dzB) &
                        &+Ppt%dx*Ppt%dy*Ppt%CvtBz%zy/(Ppt%dyF+2*Ppt%dy+Ppt%dyB)

        ! 2-O Bloque izquierdo-centro-superior
        ! It sets the column position on the matrix
        Stencil%idx_clmns(MatStencil_i,2) = i-1
        Stencil%idx_clmns(MatStencil_j,2) = j
        Stencil%idx_clmns(MatStencil_k,2) = k-1
        Stencil%Values(2)=Ppt%dy*Ppt%dz*Ppt%CvtBx%xz/(Ppt%dzF+2*Ppt%dz+Ppt%dzB) &
                        &+Ppt%dx*Ppt%dy*Ppt%CvtBz%zx/(Ppt%dxF+2*Ppt%dx+Ppt%dxB)

        ! 3-J Bloque centro-centro-superior
        ! It sets the column position on the matrix
        Stencil%idx_clmns(MatStencil_i,3) = i
        Stencil%idx_clmns(MatStencil_j,3) = j
        Stencil%idx_clmns(MatStencil_k,3) = k-1
        Stencil%Values(3)=Ppt%dy*Ppt%dz*(Ppt%CvtBx%xz-Ppt%CvtFx%xz)/(Ppt%dzF+2*Ppt%dz+Ppt%dzB) &
                        &+Ppt%dx*Ppt%dz*(Ppt%CvtBy%yz-Ppt%CvtFy%yz)/(Ppt%dzF+2*Ppt%dz+Ppt%dzB) &
                        &+Ppt%dx*Ppt%dy*2*Ppt%CvtBz%zz/(Ppt%dzB+Ppt%dz)

        ! 4-H Bloque derecho-centro-superior
        ! It sets the column position on the matrix
        Stencil%idx_clmns(MatStencil_i,4) = i+1
        Stencil%idx_clmns(MatStencil_j,4) = j
        Stencil%idx_clmns(MatStencil_k,4) = k-1
        Stencil%Values(4)=Ppt%dy*Ppt%dz*(-1)*Ppt%CvtFx%xz/(Ppt%dzF+2*Ppt%dz+Ppt%dzB) &
                        &+Ppt%dx*Ppt%dy*(-1)*Ppt%CvtBz%zx/(Ppt%dxF+2*Ppt%dx+Ppt%dxB)

        ! 5-Q Bloque centro-frontal-superior
        ! It sets the column position on the matrix
        Stencil%idx_clmns(MatStencil_i,5) = i
        Stencil%idx_clmns(MatStencil_j,5) = j+1
        Stencil%idx_clmns(MatStencil_k,5) = k-1
        Stencil%Values(5)=Ppt%dy*Ppt%dz*(-1)*Ppt%CvtFy%yz/(Ppt%dzF+2*Ppt%dz+Ppt%dzB) &
                        &+Ppt%dx*Ppt%dy*(-1)*Ppt%CvtBz%zy/(Ppt%dyF+2*Ppt%dy+Ppt%dyB)

        ! 6-M Bloque izquierdo-detras-centro
        ! It sets the column position on the matrix
        Stencil%idx_clmns(MatStencil_i,6) = i-1
        Stencil%idx_clmns(MatStencil_j,6) = j-1
        Stencil%idx_clmns(MatStencil_k,6) = k
        Stencil%Values(6)=Ppt%dy*Ppt%dz*Ppt%CvtBx%xy/(Ppt%dyF+2*Ppt%dy+Ppt%dyB) &
                        &+Ppt%dx*Ppt%dz*Ppt%CvtBy%yx/(Ppt%dxF+2*Ppt%dxB+Ppt%dx) ! subindices del diferencial malos en el paper(?)

        ! 7-F Bloque centro-detras-centro
        ! It sets the column position on the matrix
        Stencil%idx_clmns(MatStencil_i,7) = i
        Stencil%idx_clmns(MatStencil_j,7) = j-1
        Stencil%idx_clmns(MatStencil_k,7) = k
        Stencil%Values(7)=Ppt%dy*Ppt%dz* (Ppt%CvtBx%xy-Ppt%CvtFx%xy)/(Ppt%dyF+2*Ppt%dy+Ppt%dyB)  &
                        &+Ppt%dx*Ppt%dz*2*Ppt%CvtBy%yy              /(Ppt%dy+Ppt%dyB)            &
                        &+Ppt%dx*Ppt%dy* (Ppt%CvtBz%zy-Ppt%CvtFz%zy)/(Ppt%dyF+2*Ppt%dy+Ppt%dyB) 

        ! 8-D Bloque derecho-detras-centro
        ! It sets the column position on the matrix
        Stencil%idx_clmns(MatStencil_i,8) = i+1
        Stencil%idx_clmns(MatStencil_j,8) = j-1
        Stencil%idx_clmns(MatStencil_k,8) = k
        Stencil%Values(8)=Ppt%dy*Ppt%dz*(-1)*Ppt%CvtFx%xy/(Ppt%dyF+2*Ppt%dyB+Ppt%dy) & ! subindices del diferencial malos en el paper(?)
                        &+Ppt%dx*Ppt%dz*(-1)*Ppt%CvtBy%yx/(Ppt%dxF+2*Ppt%dxB+Ppt%dx)   ! subindices del diferencial malos en el paper(?)

        ! 9-K Bloque izquierdo-centro-centro
        ! It sets the column position on the matrix
        Stencil%idx_clmns(MatStencil_i,9) = i-1
        Stencil%idx_clmns(MatStencil_j,9) = j
        Stencil%idx_clmns(MatStencil_k,9) = k
        Stencil%Values(9)=Ppt%dy*Ppt%dz*2*Ppt%CvtBx%xx              /(Ppt%dx+Ppt%dxB)            &
                        &+Ppt%dx*Ppt%dz* (Ppt%CvtBy%yx-Ppt%CvtFy%yx)/(Ppt%dxF+2*Ppt%dx+Ppt%dxB)  &
                        &+Ppt%dx*Ppt%dy* (Ppt%CvtBy%zx-Ppt%CvtFy%zx)/(Ppt%dxF+2*Ppt%dx+Ppt%dxB)

        ! 10-B Bloque centro-centro-centro
        Stencil%Values(10)=Ppt%dy*Ppt%dz*(-2)*(Ppt%CvtFx%xx/(Ppt%dxF+Ppt%dx)+Ppt%CvtBx%xx/(Ppt%dx+Ppt%dxB)) &
                         &+Ppt%dx*Ppt%dz*(-2)*(Ppt%CvtFy%yy/(Ppt%dyF+Ppt%dy)+Ppt%CvtBy%yy/(Ppt%dy+Ppt%dyB)) &
                         &+Ppt%dx*Ppt%dy*(-2)*(Ppt%CvtFz%zz/(Ppt%dzF+Ppt%dz)+Ppt%CvtBz%zz/(Ppt%dz+Ppt%dzB))

        ! 11-A Bloque derecho-centro-centro
        ! It sets the column position on the matrix
        Stencil%idx_clmns(MatStencil_i,11) = i+1
        Stencil%idx_clmns(MatStencil_j,11) = j
        Stencil%idx_clmns(MatStencil_k,11) = k
        Stencil%Values(11)=Ppt%dy*Ppt%dz*2*Ppt%CvtFx%xx             /(Ppt%dxF+Ppt%dx)            &
                        &+Ppt%dx*Ppt%dz* (Ppt%CvtFy%yx-Ppt%CvtBy%yx)/(Ppt%dxF+2*Ppt%dx+Ppt%dxB)  &
                        &+Ppt%dx*Ppt%dy* (Ppt%CvtFz%zx-Ppt%CvtBz%zx)/(Ppt%dxF+2*Ppt%dx+Ppt%dxB)

        ! 12-L Bloque izquierdo-frontal-centro
        ! It sets the column position on the matrix
        Stencil%idx_clmns(MatStencil_i,12) = i-1
        Stencil%idx_clmns(MatStencil_j,12) = j+1
        Stencil%idx_clmns(MatStencil_k,12) = k
        Stencil%Values(12)=Ppt%dy*Ppt%dz*(-1)*Ppt%CvtBx%xy/(Ppt%dyF+2*Ppt%dy+Ppt%dyB) &
                         &+Ppt%dx*Ppt%dz*(-1)*Ppt%CvtFy%yx/(Ppt%dxF+2*Ppt%dxB+Ppt%dx)   ! subindices del diferencial malos en el paper(?)

        ! 13-E Bloque centro-frontal-centro
        ! It sets the column position on the matrix
        Stencil%idx_clmns(MatStencil_i,13) = i
        Stencil%idx_clmns(MatStencil_j,13) = j+1
        Stencil%idx_clmns(MatStencil_k,13) = k
        Stencil%Values(13)=Ppt%dy*Ppt%dz* (Ppt%CvtFx%xy-Ppt%CvtBx%xy)/(Ppt%dyF+2*Ppt%dy+Ppt%dyB)  &
                         &+Ppt%dx*Ppt%dz*2*Ppt%CvtFy%yy              /(Ppt%dyF+Ppt%dy)            &
                         &+Ppt%dx*Ppt%dy* (Ppt%CvtFz%zy-Ppt%CvtBz%zy)/(Ppt%dyF+2*Ppt%dy+Ppt%dyB) 

        ! 14-C Bloque derecho-frontal-centro
        ! It sets the column position on the matrix
        Stencil%idx_clmns(MatStencil_i,14) = i+1
        Stencil%idx_clmns(MatStencil_j,14) = j+1
        Stencil%idx_clmns(MatStencil_k,14) = k
        Stencil%Values(14)=Ppt%dy*Ppt%dz*Ppt%CvtFx%xy/(Ppt%dyF+2*Ppt%dyB+Ppt%dy) & ! subindices del diferencial malos en el paper(?)
                         &+Ppt%dx*Ppt%dz*Ppt%CvtFy%yx/(Ppt%dxF+2*Ppt%dxB+Ppt%dx)   ! subindices del diferencial malos en el paper(?)

        ! 15-R Bloque centro-detras-inferior
        ! It sets the column position on the matrix
        Stencil%idx_clmns(MatStencil_i,15) = i
        Stencil%idx_clmns(MatStencil_j,15) = j-1
        Stencil%idx_clmns(MatStencil_k,15) = k+1
        Stencil%Values(15)=Ppt%dy*Ppt%dz*(-1)*Ppt%CvtBy%yz/(Ppt%dzF+2*Ppt%dz+Ppt%dzB) &
                         &+Ppt%dx*Ppt%dy*(-1)*Ppt%CvtFz%zy/(Ppt%dyF+2*Ppt%dy+Ppt%dyB)

        ! 16-N Bloque izquierdo-centro-inferior
        ! It sets the column position on the matrix
        Stencil%idx_clmns(MatStencil_i,16) = i-1
        Stencil%idx_clmns(MatStencil_j,16) = j
        Stencil%idx_clmns(MatStencil_k,16) = k+1
        Stencil%Values(16)=Ppt%dy*Ppt%dz*(-1)*Ppt%CvtBx%xz/(Ppt%dzF+2*Ppt%dz+Ppt%dzB) &
                         &+Ppt%dx*Ppt%dy*(-1)*Ppt%CvtFz%zx/(Ppt%dxF+2*Ppt%dx+Ppt%dxB)

        ! 17-I Bloque centro-centro-inferior
        ! It sets the column position on the matrix
        Stencil%idx_clmns(MatStencil_i,17) = i
        Stencil%idx_clmns(MatStencil_j,17) = j
        Stencil%idx_clmns(MatStencil_k,17) = k+1
        Stencil%Values(17)=Ppt%dy*Ppt%dz*(Ppt%CvtFx%xz-Ppt%CvtBx%xz)/(Ppt%dzF+2*Ppt%dz+Ppt%dzB) &
                         &+Ppt%dx*Ppt%dz*(Ppt%CvtFy%yz-Ppt%CvtBy%yz)/(Ppt%dzF+2*Ppt%dz+Ppt%dzB) &
                         &+Ppt%dx*Ppt%dy*2*Ppt%CvtFz%zz/(Ppt%dzF+Ppt%dz)

        ! 18-G Bloque derecho-centro-inferior
        ! It sets the column position on the matrix
        Stencil%idx_clmns(MatStencil_i,18) = i+1
        Stencil%idx_clmns(MatStencil_j,18) = j
        Stencil%idx_clmns(MatStencil_k,18) = k+1
        Stencil%Values(18)=Ppt%dy*Ppt%dz*Ppt%CvtFx%xz/(Ppt%dzF+2*Ppt%dz+Ppt%dzB) &
                         &+Ppt%dx*Ppt%dy*Ppt%CvtFz%zx/(Ppt%dxF+2*Ppt%dx+Ppt%dxB)

        ! 19-P Bloque centro-frontal-inferior
        ! It sets the column position on the matrix
        Stencil%idx_clmns(MatStencil_i,19) = i
        Stencil%idx_clmns(MatStencil_j,19) = j+1
        Stencil%idx_clmns(MatStencil_k,19) = k+1
        Stencil%Values(19)=Ppt%dy*Ppt%dz*Ppt%CvtFy%yz/(Ppt%dzF+2*Ppt%dz+Ppt%dzB) &
                         &+Ppt%dx*Ppt%dy*Ppt%CvtFz%zy/(Ppt%dyF+2*Ppt%dy+Ppt%dyB)

    ELSEIF ((Ppt%StnclTplgy(10).EQ.2).OR.(Ppt%StnclTplgy(10).EQ.0)) THEN ! Dirichlet or inactive cell
        Stencil%Values(10)=-one
    END IF

END SUBROUTINE GetLiStencil

! SUBROUTINE AnisoflowStencil(Ppt,Stencil,ierr)

!     USE ANISOFLOW_Types, ONLY : Property,StencilVar

!     IMPLICIT NONE

! #include <petsc/finclude/petscsys.h>

!     PetscErrorCode,INTENT(INOUT)            :: ierr
!     TYPE(Property),INTENT(IN)               :: Ppt
!     TYPE(StencilVar),INTENT(INOUT)          :: Stencil

!     PetscInt                                :: i,j,k
!     PetscReal                               :: one=1.d0

!     ! The Li stencil is based on 19 cells.
!     ALLOCATE(Stencil%idx_clmns(4,19))
!     ALLOCATE(Stencil%Values(19))
!     Stencil%Size=19

!     ! It gets the position of the cell.
!     i=Ppt%Pstn%i
!     j=Ppt%Pstn%j
!     k=Ppt%Pstn%k

!     ! Rows to modify
!     Stencil%idx_rws(MatStencil_i,1) = i
!     Stencil%idx_rws(MatStencil_j,1) = j
!     Stencil%idx_rws(MatStencil_k,1) = k
!     ! Columns to modify
!     Stencil%idx_clmns(MatStencil_i,:) = i
!     Stencil%idx_clmns(MatStencil_j,:) = j
!     Stencil%idx_clmns(MatStencil_k,:) = k

!     ! Initial stencil vaulues
!     Stencil%Values(:)=0.0

!     ! If the current cell is an active cell:
!     IF (Ppt%StnclTplgy(10).EQ.1) THEN

!         ! 1-S Bloque centro-detras-superior
!         ! It sets the column position on the matrix
!         Stencil%idx_clmns(MatStencil_i,1) = i
!         Stencil%idx_clmns(MatStencil_j,1) = j-1
!         Stencil%idx_clmns(MatStencil_k,1) = k-1
!         Stencil%Values(1)=

!         ! 2-O Bloque izquierdo-centro-superior
!         ! It sets the column position on the matrix
!         Stencil%idx_clmns(MatStencil_i,2) = i-1
!         Stencil%idx_clmns(MatStencil_j,2) = j
!         Stencil%idx_clmns(MatStencil_k,2) = k-1
!         Stencil%Values(2)=

!         ! 3-J Bloque centro-centro-superior
!         ! It sets the column position on the matrix
!         Stencil%idx_clmns(MatStencil_i,3) = i
!         Stencil%idx_clmns(MatStencil_j,3) = j
!         Stencil%idx_clmns(MatStencil_k,3) = k-1
!         Stencil%Values(3)=Ppt%dy*Ppt%dz*KxxB/Ppt%dxB 

!         ! 4-H Bloque derecho-centro-superior
!         ! It sets the column position on the matrix
!         Stencil%idx_clmns(MatStencil_i,4) = i+1
!         Stencil%idx_clmns(MatStencil_j,4) = j
!         Stencil%idx_clmns(MatStencil_k,4) = k-1
!         Stencil%Values(4)=

!         ! 5-Q Bloque centro-frontal-superior
!         ! It sets the column position on the matrix
!         Stencil%idx_clmns(MatStencil_i,5) = i
!         Stencil%idx_clmns(MatStencil_j,5) = j+1
!         Stencil%idx_clmns(MatStencil_k,5) = k-1
!         Stencil%Values(5)=

!         ! 6-M Bloque izquierdo-detras-centro
!         ! It sets the column position on the matrix
!         Stencil%idx_clmns(MatStencil_i,6) = i-1
!         Stencil%idx_clmns(MatStencil_j,6) = j-1
!         Stencil%idx_clmns(MatStencil_k,6) = k
!         Stencil%Values(6)=-(Ppt%dy*Ppt%dz(CvtBx%xy/()))

!         ! 7-F Bloque centro-detras-centro
!         ! It sets the column position on the matrix
!         Stencil%idx_clmns(MatStencil_i,7) = i
!         Stencil%idx_clmns(MatStencil_j,7) = j-1
!         Stencil%idx_clmns(MatStencil_k,7) = k
!         Stencil%Values(7)=Ppt%dx*Ppt%dz*KyyB/Ppt%dyB 

!         ! 8-D Bloque derecho-detras-centro
!         ! It sets the column position on the matrix
!         Stencil%idx_clmns(MatStencil_i,8) = i+1
!         Stencil%idx_clmns(MatStencil_j,8) = j-1
!         Stencil%idx_clmns(MatStencil_k,8) = k
!         Stencil%Values(8)=

!         ! 9-K Bloque izquierdo-centro-centro
!         ! It sets the column position on the matrix
!         Stencil%idx_clmns(MatStencil_i,9) = i-1
!         Stencil%idx_clmns(MatStencil_j,9) = j
!         Stencil%idx_clmns(MatStencil_k,9) = k
!         Stencil%Values(9)=Ppt%dx*Ppt%dy*KzzB/Ppt%dzB

!         ! 10-B Bloque centro-centro-centro
!         Stencil%Values(10)=-(Ppt%dy*Ppt%dz*KxxF/Ppt%dxF+Ppt%dy*Ppt%dz*KxxB/Ppt%dxB &
!                         & +Ppt%dx*Ppt%dz*KyyF/Ppt%dyF+Ppt%dx*Ppt%dz*KyyB/Ppt%dyB              &
!                         & +Ppt%dx*Ppt%dy*KzzF/Ppt%dzF+Ppt%dx*Ppt%dy*KzzB/Ppt%dzB              &
!                         & +Ss*Ppt%dx*Ppt%dy*Ppt%dz)

!         ! 11-A Bloque derecho-centro-centro
!         ! It sets the column position on the matrix
!         Stencil%idx_clmns(MatStencil_i,11) = i+1
!         Stencil%idx_clmns(MatStencil_j,11) = j
!         Stencil%idx_clmns(MatStencil_k,11) = k
!         Stencil%Values(11)=pt%dy*Ppt%dz*KxxF/Ppt%dxF

!         ! 12-L Bloque izquierdo-frontal-centro
!         ! It sets the column position on the matrix
!         Stencil%idx_clmns(MatStencil_i,12) = i-1
!         Stencil%idx_clmns(MatStencil_j,12) = j+1
!         Stencil%idx_clmns(MatStencil_k,12) = k
!         Stencil%Values(12)=

!         ! 13-E Bloque centro-frontal-centro
!         ! It sets the column position on the matrix
!         Stencil%idx_clmns(MatStencil_i,13) = i
!         Stencil%idx_clmns(MatStencil_j,13) = j+1
!         Stencil%idx_clmns(MatStencil_k,13) = k
!         Stencil%Values(13)=Ppt%dx*Ppt%dz*KyyF/Ppt%dyF

!         ! 14-C Bloque derecho-frontal-centro
!         ! It sets the column position on the matrix
!         Stencil%idx_clmns(MatStencil_i,14) = i+1
!         Stencil%idx_clmns(MatStencil_j,14) = j+1
!         Stencil%idx_clmns(MatStencil_k,14) = k
!         Stencil%Values(14)=

!         ! 15-R Bloque centro-detras-inferior
!         ! It sets the column position on the matrix
!         Stencil%idx_clmns(MatStencil_i,15) = i
!         Stencil%idx_clmns(MatStencil_j,15) = j-1
!         Stencil%idx_clmns(MatStencil_k,15) = k+1
!         Stencil%Values(15)=

!         ! 16-N Bloque izquierdo-centro-inferior
!         ! It sets the column position on the matrix
!         Stencil%idx_clmns(MatStencil_i,16) = i-1
!         Stencil%idx_clmns(MatStencil_j,16) = j
!         Stencil%idx_clmns(MatStencil_k,16) = k+1
!         Stencil%Values(16)=

!         ! 17-I Bloque centro-centro-inferior
!         ! It sets the column position on the matrix
!         Stencil%idx_clmns(MatStencil_i,17) = i
!         Stencil%idx_clmns(MatStencil_j,17) = j
!         Stencil%idx_clmns(MatStencil_k,17) = k+1
!         Stencil%Values(17)=Ppt%dx*Ppt%dy*KzzF/Ppt%dzF

!         ! 18-G Bloque derecho-centro-inferior
!         ! It sets the column position on the matrix
!         Stencil%idx_clmns(MatStencil_i,18) = i+1
!         Stencil%idx_clmns(MatStencil_j,18) = j
!         Stencil%idx_clmns(MatStencil_k,18) = k+1
!         Stencil%Values(18)=

!         ! 19-P Bloque centro-frontal-inferior
!         ! It sets the column position on the matrix
!         Stencil%idx_clmns(MatStencil_i,19) = i
!         Stencil%idx_clmns(MatStencil_j,19) = j+1
!         Stencil%idx_clmns(MatStencil_k,19) = k+1
!         Stencil%Values(19)=        

!     ! If the current cell is an Neumman x cell:
!     ELSEIF (Ppt%StnclTplgy(10).EQ.3) THEN ! Neumman x

!         Stencil%Values(10)=-one

!         ! It is a decision of wich cell is in Neumman x condition to assign the equality.
!         IF ((Ppt%StnclTplgy(9).EQ.1).AND.(Ppt%StnclTplgy(11).EQ.1)) THEN
!             CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
!                 & "ERROR: Neumman on x axis was bad defined\n",ierr)
!             STOP
!         ELSEIF (Ppt%StnclTplgy(9).EQ.1) THEN
!             Stencil%Values(9)=one
!             Stencil%idx_clmns(MatStencil_i,9) = i-1
!             Stencil%idx_clmns(MatStencil_j,9) = j
!             Stencil%idx_clmns(MatStencil_k,9) = k
!         ELSEIF (Ppt%StnclTplgy(11).EQ.1) THEN
!             Stencil%Values(11)=one
!             Stencil%idx_clmns(MatStencil_i,11) = i+1
!             Stencil%idx_clmns(MatStencil_j,11) = j
!             Stencil%idx_clmns(MatStencil_k,11) = k
!         END IF

!     ! If the current cell is an Neumman y cell:
!     ELSEIF (Ppt%StnclTplgy(10).EQ.4) THEN ! Neumman y
!         Stencil%Values(10)=-one
!         ! It is a decision of wich cell is in Neumman y condition to assign the equality.
!         IF ((Ppt%StnclTplgy(7).EQ.1).AND.(Ppt%StnclTplgy(13).EQ.1)) THEN
!             CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
!                 & "ERROR: Neumman on y axis was bad defined\n",ierr)
!             STOP
!         ELSEIF (Ppt%StnclTplgy(7).EQ.1) THEN
!             Stencil%Values(7)=one
!             Stencil%idx_clmns(MatStencil_i,7) = i
!             Stencil%idx_clmns(MatStencil_j,7) = j-1
!             Stencil%idx_clmns(MatStencil_k,7) = k
!         ELSEIF (Ppt%StnclTplgy(13).EQ.1) THEN
!             Stencil%Values(13)=one
!             Stencil%idx_clmns(MatStencil_i,13) = i
!             Stencil%idx_clmns(MatStencil_j,13) = j+1
!             Stencil%idx_clmns(MatStencil_k,13) = k
!         END IF

!     ! If the current cell is an Neumman z cell:
!     ELSEIF (Ppt%StnclTplgy(10).EQ.5) THEN ! Neumman z
!         Stencil%Values(10)=-one
!         ! It is a decision of wich cell is in Neumman z condition to assign the equality.
!         IF ((Ppt%StnclTplgy(3).EQ.1).AND.(Ppt%StnclTplgy(17).EQ.1)) THEN
!             CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,                         &
!                 & "ERROR: Neumman on z axis was bad defined\n",ierr)
!             STOP
!         ELSEIF (Ppt%StnclTplgy(3).EQ.1) THEN
!             Stencil%Values(3)=one
!             Stencil%idx_clmns(MatStencil_i,3) = i
!             Stencil%idx_clmns(MatStencil_j,3) = j
!             Stencil%idx_clmns(MatStencil_k,3) = k-1
!         ELSEIF (Ppt%StnclTplgy(17).EQ.1) THEN
!             Stencil%Values(17)=one
!             Stencil%idx_clmns(MatStencil_i,17) = i
!             Stencil%idx_clmns(MatStencil_j,17) = j
!             Stencil%idx_clmns(MatStencil_k,17) = k+1
!         END IF
!     ELSEIF (Ppt%StnclTplgy(10).EQ.2) THEN ! Dirichlet
!         Stencil%Values(10)=-one
!     END IF

! END SUBROUTINE AnisoflowStencil

SUBROUTINE ApplyDirichlet(BCFld,TimeZone,b,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,BoundaryConditions
    USe ANISOFLOW_Interface, ONLY : GetVerbose

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(BoundaryConditions),INTENT(IN)     :: BCFld
    PetscInt,INTENT(IN)                     :: TimeZone
    Vec,INTENT(INOUT)                       :: b

    VecScatter                              :: Scatter
    IS                                      :: NaturalOrder
    PetscInt                                :: DirichSize
    PetscBool                               :: Verbose

    CALL VecGetSize(BCFld%Dirich(TimeZone),DirichSize,ierr)

    CALL ISCreateStride(PETSC_COMM_WORLD,DirichSize,0,1,NaturalOrder,ierr)
    CALL VecScatterCreate(BCFld%Dirich(TimeZone),NaturalOrder,b,BCFld%DirichIS(TimeZone),Scatter,ierr)

    CALL VecScatterBegin(Scatter,BCFld%Dirich(TimeZone),b,INSERT_VALUES,SCATTER_FORWARD,ierr)
    CALL VecScatterEnd(Scatter,BCFld%Dirich(TimeZone),b,INSERT_VALUES,SCATTER_FORWARD,ierr)

    CALL VecScatterDestroy(Scatter,ierr)
    CALL ISDestroy(NaturalOrder,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetSystem Event] Dirichlet boundary conditions properly implemented\n",ierr)

END SUBROUTINE ApplyDirichlet

SUBROUTINE ApplySource(BCFld,TimeZone,b,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,BoundaryConditions
    USe ANISOFLOW_Interface, ONLY : GetVerbose

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(BoundaryConditions),INTENT(IN)     :: BCFld
    PetscInt,INTENT(IN)                     :: TimeZone
    Vec,INTENT(INOUT)                       :: b

    VecScatter                              :: Scatter
    IS                                      :: NaturalOrder
    PetscInt                                :: SourceSize
    PetscBool                               :: Verbose

    CALL VecGetSize(BCFld%Source(TimeZone),SourceSize,ierr)

    CALL ISCreateStride(PETSC_COMM_WORLD,SourceSize,0,1,NaturalOrder,ierr)
    CALL VecScatterCreate(BCFld%Source(TimeZone),NaturalOrder,b,BCFld%SourceIS(TimeZone),Scatter,ierr)

    CALL VecScatterBegin(Scatter,BCFld%Source(TimeZone),b,ADD_VALUES,SCATTER_FORWARD,ierr)
    CALL VecScatterEnd(Scatter,BCFld%Source(TimeZone),b,ADD_VALUES,SCATTER_FORWARD,ierr)

    CALL VecScatterDestroy(Scatter,ierr)
    CALL ISDestroy(NaturalOrder,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetSystem Event] Source terms properly implemented\n",ierr)

END SUBROUTINE ApplySource

SUBROUTINE ApplyCauchy(BCFld,TimeZone,A,b,ierr)

    USE ANISOFLOW_Types, ONLY : Geometry,BoundaryConditions
    USe ANISOFLOW_Interface, ONLY : GetVerbose

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(BoundaryConditions),INTENT(IN)     :: BCFld
    PetscInt,INTENT(IN)                     :: TimeZone
    Vec,INTENT(INOUT)                       :: b
    Mat,INTENT(INOUT)                       :: A

    Vec                                     :: MultHC,MultHeC,DiagA,DiagApartial
    VecScatter                              :: Scatter
    IS                                      :: NaturalOrder
    PetscBool                               :: Verbose

!     CALL VecView(BCFld%CauchyHe(TimeZone),PETSC_VIEWER_STDOUT_WORLD,ierr)
!     CALL VecView(BCFld%CauchyC(TimeZone),PETSC_VIEWER_STDOUT_WORLD,ierr)

    CALL VecDuplicate(BCFld%CauchyHe(TimeZone),MultHC,ierr)
    CALL VecDuplicate(BCFld%CauchyHe(TimeZone),MultHeC,ierr)
    CALL VecDuplicate(BCFld%CauchyHe(TimeZone),DiagApartial,ierr)

    CALL VecPointwiseMult(MultHeC,BCFld%CauchyC(TimeZone),BCFld%CauchyHe(TimeZone),ierr)

!     CALL VecView(MultHeC,PETSC_VIEWER_STDOUT_WORLD,ierr)

    CALL ISCreateStride(PETSC_COMM_WORLD,BCFld%SizeCauchy,0,1,NaturalOrder,ierr)
    CALL VecScatterCreate(MultHeC,NaturalOrder,b,BCFld%CauchyIS(TimeZone),Scatter,ierr)

    CALL VecScatterBegin(Scatter,MultHeC,b,ADD_VALUES,SCATTER_FORWARD,ierr)
    CALL VecScatterEnd(Scatter,MultHeC,b,ADD_VALUES,SCATTER_FORWARD,ierr)

    CALL VecDestroy(MultHeC,ierr)

    CALL VecDuplicate(b,DiagA,ierr)
    CALL MatGetDiagonal(A,DiagA,ierr)

    CALL VecScatterBegin(Scatter,DiagA,DiagApartial,ADD_VALUES,SCATTER_REVERSE,ierr)
    CALL VecScatterEnd(Scatter,DiagA,DiagApartial,ADD_VALUES,SCATTER_REVERSE,ierr)

!     CALL VecView(DiagA,PETSC_VIEWER_STDOUT_WORLD,ierr)
!     CALL VecView(DiagApartial,PETSC_VIEWER_STDOUT_WORLD,ierr)

    CALL VecPointwiseMult(DiagApartial,BCFld%CauchyC(TimeZone),MultHC,ierr)

    CALL VecScatterBegin(Scatter,DiagApartial,DiagA,ADD_VALUES,SCATTER_FORWARD,ierr)
    CALL VecScatterEnd(Scatter,DiagApartial,DiagA,ADD_VALUES,SCATTER_FORWARD,ierr)

    CALL MatDiagonalSet(A,diagA,INSERT_VALUES,ierr)

    CALL VecScatterDestroy(Scatter,ierr)
    CALL ISDestroy(NaturalOrder,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetSystem Event] Cauchy terms properly implemented\n",ierr)

END SUBROUTINE ApplyCauchy

SUBROUTINE ApplyTimeDiff(PptFld,BCFld,TimeZone,TimeStep,A,b,x,ierr)

    USE ANISOFLOW_Types, ONLY : PropertiesField,BoundaryConditions
    USE ANISOFLOW_Interface, ONLY : GetVerbose

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(PropertiesField),INTENT(IN)        :: PptFld
    TYPE(BoundaryConditions),INTENT(IN)     :: BCFld
    PetscInt,INTENT(IN)                     :: TimeZone,TimeStep
    Mat,INTENT(INOUT)                       :: A
    Vec,INTENT(INOUT)                       :: b,x

    Vec                                     :: VecDT
    PetscReal                               :: DT,one=1.d0
    CHARACTER(LEN=200)                      :: CharDT
    PetscBool                               :: Verbose

    CALL GetVerbose(Verbose,ierr)
    CALL GetDT(BCFld,TimeZone,TimeStep,DT,ierr)
    WRITE(CharDT,*)DT
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[GetSystem Stage] DT: "//ADJUSTL(TRIM(CharDT))//"\n",ierr)

    CALL MatGetDiagonal(A,VecDT,ierr)
    CALL VecView(VecDT,PETSC_VIEWER_STDOUT_WORLD,ierr)
    STOP

    CALL VecDuplicate(x,VecDT,ierr)
    CALL VecSet(VecDT,-one/DT,ierr)
    CALL VecPointwiseMult(VecDT,VecDT,PptFld%Sto%Cell,ierr)
    CALL MatDiagonalSet(A,VecDT,ADD_VALUES,ierr)
    CALL VecPointwiseMult(b,VecDT,x,ierr)
    CALL VecAXPY(b,one,VecDT,ierr)
    CALL VecDestroy(VecDT,ierr)

!     CALL VecDuplicate(x,VecDT,ierr)
!     CALL VecSet(VecDT,-one/DT,ierr)
!     CALL MatDiagonalSet(A,VecDT,ADD_VALUES,ierr)
!     CALL VecDestroy(VecDT,ierr)

!     CALL VecAXPY(b,-one/DT,x,ierr)



END SUBROUTINE ApplyTimeDiff

SUBROUTINE GetDT(BCFld,TimeZone,TimeStep,DT,ierr)

    USE ANISOFLOW_Types, ONLY : BoundaryConditions

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(BoundaryConditions),INTENT(IN)     :: BCFld
    PetscInt,INTENT(IN)                     :: TimeZone,TimeStep
    PetscReal,INTENT(OUT)                   :: DT

    PetscInt                                :: TimeZoneSize
    PetscReal                               :: Time1,Time2
    !PetscReal,POINTER                       :: TimeArray(:)

    Time2=BCFld%TimeZone(TimeZone)%Time(TimeStep)
    IF ((TimeStep.EQ.1).AND.(TimeZone.EQ.1)) THEN
        Time1=0.d0
    ELSEIF (TimeStep.EQ.1) THEN
        TimeZoneSize=SIZE(BCFld%TimeZone(TimeZone-1)%Time) 
        Time1=BCFld%TimeZone(TimeZone-1)%Time(TimeZoneSize)
        Time1=0.d0
    ELSE 
        Time1=BCFld%TimeZone(TimeZone)%Time(TimeStep-1)
    END IF


    DT=Time2-Time1

END SUBROUTINE GetDT

SUBROUTINE DestroySystem(A,b,x,ierr)

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>

    PetscErrorCode,INTENT(INOUT)        :: ierr
    Mat,INTENT(INOUT)                   :: A
    Vec,INTENT(INOUT)                   :: b,x
    CALL VecDestroy(b,ierr)
    CALL VecDestroy(x,ierr)
    CALL MatDestroy(A,ierr)

END SUBROUTINE DestroySystem

END MODULE ANISOFLOW_BuildSystem