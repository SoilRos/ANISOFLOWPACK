MODULE ANISOFLOW_BuildSystem

    IMPLICIT NONE

CONTAINS

SUBROUTINE BuildSystem(Gmtry,PptFld,A,ierr)

    USE ANISOFLOW_Types,      ONLY : Geometry,PropertiesField,       &
                                   & Property,StencilVar
    USE ANISOFLOW_Properties, ONLY : GetLocalProperty
    USE ANISOFLOW_Interface,  ONLY : GetVerbose

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
    CALL PetscLogEventBegin(Event,PETSC_NULL_OBJECT,                 &
        & PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,      &
        & "["//ADJUSTL(TRIM(EventName))//" Event] Inizialited\n",ierr)
    

    CALL DMCreateMatrix(Gmtry%DataMngr,A,ierr)

    CALL DMDAGetCorners(Gmtry%DataMngr,corn(1),corn(2),corn(3),      &
        & widthL(1),widthL(2),widthL(3),ierr)

    DO k=corn(3),corn(3)+widthL(3)-1
        DO j=corn(2),corn(2)+widthL(2)-1
            DO i=corn(1),corn(1)+widthL(1)-1

                CALL GetLocalProperty(Gmtry,PptFld,Ppt,i,j,k,ierr) 
                CALL GetStencil(Ppt,Stencil,ierr)
                
                IF (Stencil%IsActive) THEN
                    CALL MatSetValuesStencil(A,1,Stencil%idx_rws,    &
                        & Stencil%idx_size,Stencil%idx_clmns,        &
                        & Stencil%idx_val,ADD_VALUES,ierr)
                END IF
            END DO
        END DO
    END DO

    CALL MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)  

    CALL MatZeroRowsColumnsIS(A,Gmtry%InactiveIS,-1.D0,              &
        & PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,      &
        & "["//ADJUSTL(TRIM(EventName))//" Event] Finalized\n",ierr)
    
    CALL PetscLogFlops(EventFlops,ierr)
    CALL PetscLogEventEnd(Event,PETSC_NULL_OBJECT,                   &
        & PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

END SUBROUTINE BuildSystem

SUBROUTINE GetInitSol(Gmtry,x,ierr)

    USE ANISOFLOW_Types,     ONLY : InputTypeVar,Geometry
    USE ANISOFLOW_Interface, ONLY : GetInputDir,GetInputTypeInitSol, &
                                  & GetInputFileInitSol,             &
                                  & GetInitSolUniValue,GetVerbose

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscviewer.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Geometry),INTENT(IN)               :: Gmtry
    Vec,INTENT(OUT)                         :: x

    PetscBool                   :: InputFileInitSolFlg,              &
                                 & InitSolUniValueFlg,Verbose
    PetscReal                   :: InitSolUniValue
    TYPE(InputTypeVar)          :: InputType
    CHARACTER(LEN=200)          :: Route,InputDir,InputFileInitSol,  &
                                 & CharInitSolUniValue
    PetscViewer                 :: Viewer
   
    CALL GetVerbose(Verbose,ierr)
    CALL GetInputDir(InputDir,ierr)
    CALL GetInputFileInitSol(InputFileInitSol,InputFileInitSolFlg,   &
        & ierr)

    CALL DMCreateGlobalVector(Gmtry%DataMngr,x,ierr)

    IF (InputFileInitSolFlg) THEN
        CALL GetInputTypeInitSol(InputType,ierr)
        Route=ADJUSTL(TRIM(InputDir)//TRIM(InputFileInitSol))
        CALL PetscObjectSetName(x,"Solution",ierr)
        IF (InputType%InitSol.EQ.1) THEN
            CALL PetscViewerBinaryOpen(PETSC_COMM_WORLD,Route,       &
                & FILE_MODE_READ,Viewer,ierr)
            CALL VecLoad(x,Viewer,ierr)
            CALL PetscViewerDestroy(Viewer,ierr)
            IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,&
                & "[SolveSystem Event] Inital solution "//           &
                & ADJUSTL(TRIM(InputFileInitSol))//                  &
                & " has been implemented properly.\n",ierr)
        ELSEIF (InputType%InitSol.EQ.2) THEN
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,           &
                & "ERROR: Initial solution is not able to be opend"//&
                & " from ASCII format, you should use binary or HD"//&
                & "F5 formats. Initial solution was set to 0.\n",ierr)
        ELSEIF (InputType%InitSol.EQ.3) THEN
#if defined(PETSC_HAVE_HDF5)
            CALL PetscViewerHDF5Open(PETSC_COMM_WORLD,Route,         &
                & FILE_MODE_READ,Viewer,ierr)
            CALL VecLoad(x,Viewer,ierr)
            CALL PetscViewerDestroy(Viewer,ierr)
            IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,&
                & "[SolveSystem Event] Inital solution "//           &
                & ADJUSTL(TRIM(InputFileInitSol))//                  &
                & " has been implemented properly.\n",ierr)
#else
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,           &
                & "ERROR: Initial solution is not able to be opend"//&
                & " from HDF5 format, you should use binary or ins"//&
                & "tall HDF5 libraries. Initial solution was set t"//&
                & "o 0.\n",ierr)
#endif
        ELSE 
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,           &
                & "ERROR: InitSol is not valid. Initial solution w"//&
                & "as set to 0.\n",ierr)
        END IF
    ELSE
        CALL GetInitSolUniValue(InitSolUniValue,InitSolUniValueFlg,  &
            & ierr)
        IF (InitSolUniValueFlg) THEN 
            CALL VecSet(x,InitSolUniValue,ierr)
            WRITE(CharInitSolUniValue,*)InitSolUniValue
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,           &
                & "[SolveSystem Event] Inital solution was set "//   &
                & ADJUSTL(TRIM(CharInitSolUniValue))//" .\n",ierr)
        ELSE
            CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,           &
                & "WARING: Initial solution was set to 0.\n",ierr)
        END IF
    END IF

END SUBROUTINE GetInitSol

SUBROUTINE GetStencil(Ppt,Stencil,ierr)

    USE ANISOFLOW_Types,     ONLY : Property,StencilVar,RunOptionsVar
    USE ANISOFLOW_Interface, ONLY : GetRunOptions

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Property),INTENT(IN)               :: Ppt
    TYPE(StencilVar),INTENT(OUT)            :: Stencil
    
    TYPE(RunOptionsVar)                     :: RunOptions
        
    ! Quitar esto, toma mucho tiempo en cada llamada!
    CALL GetRunOptions(RunOptions,ierr) 

    ALLOCATE(Stencil%idx_rws(4,1))

    IF (RunOptions%Scheme.EQ.1) THEN
        CALL GetTraditionalStencil(Ppt,Stencil,ierr)
    ELSEIF (RunOptions%Scheme.EQ.2) THEN
        CALL GetLiStencil(Ppt,Stencil,ierr)
    ELSEIF (RunOptions%Scheme.EQ.3) THEN
        CALL GetPerezStencil(Ppt,Stencil,ierr)
    ELSE
    
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,               &
            & "ERROR: Run_options_scheme command must be an intege"//&
            & "r between 1 and 2\n",ierr)
        STOP
    
    END IF

END SUBROUTINE GetStencil

SUBROUTINE GetTraditionalStencil(Ppt,Stencil,ierr)

    USE ANISOFLOW_Types,    ONLY : Property,StencilVar
    USE ANISOFLOW_Operators

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Property),INTENT(IN)               :: Ppt
    TYPE(StencilVar),INTENT(INOUT)          :: Stencil

    PetscInt                                :: i,j,k
    PetscReal                               :: one=1.d0,zero=0.d0

    IF (Ppt%StnclType.NE.1) THEN
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,               &
            & "ERROR: Properties must have a a stencil type star.\n",&
            & ierr)
        STOP
    END IF

    ALLOCATE(Stencil%idx_clmns(4,7))
    ALLOCATE(Stencil%idx_val(7))
    Stencil%idx_size=7

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
    Stencil%idx_val(:)=zero

    ! If the current cell is an active, source, or cauchy cell:
    IF ((Ppt%StnclTplgy(4).EQ.1).OR.(Ppt%StnclTplgy(4).EQ.3).OR.     &
        & (Ppt%StnclTplgy(4).EQ.4)) THEN

        IF (Ppt%StnclTplgy(1).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,1) = i
            Stencil%idx_clmns(MatStencil_j,1) = j
            Stencil%idx_clmns(MatStencil_k,1) = k-1
            Stencil%idx_val(1)=(Ppt%dx(4)*Ppt%dy(4))*                &
                & (Ppt%Cvt(1)%zz.ARMONIC.Ppt%Cvt(4)%zz)*             &
                & 2/(Ppt%dz(1)+Ppt%dz(4))
        END IF

        IF (Ppt%StnclTplgy(2).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,2) = i
            Stencil%idx_clmns(MatStencil_j,2) = j-1
            Stencil%idx_clmns(MatStencil_k,2) = k
            Stencil%idx_val(2)=(Ppt%dx(4)*Ppt%dz(4))*                &
                & (Ppt%Cvt(2)%yy.ARMONIC.Ppt%Cvt(4)%yy)*             &
                & 2/(Ppt%dy(2)+Ppt%dy(4))
        END IF

        IF (Ppt%StnclTplgy(3).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,3) = i-1
            Stencil%idx_clmns(MatStencil_j,3) = j
            Stencil%idx_clmns(MatStencil_k,3) = k
            Stencil%idx_val(3)=(Ppt%dy(4)*Ppt%dz(4))*                &
                & (Ppt%Cvt(3)%xx.ARMONIC.Ppt%Cvt(4)%xx)*             &
                & 2/(Ppt%dx(3)+Ppt%dx(4))
        END IF

        IF (Ppt%StnclTplgy(5).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,5) = i+1
            Stencil%idx_clmns(MatStencil_j,5) = j
            Stencil%idx_clmns(MatStencil_k,5) = k
            Stencil%idx_val(5)=(Ppt%dy(4)*Ppt%dz(4))*                &
                & (Ppt%Cvt(5)%xx.ARMONIC.Ppt%Cvt(4)%xx)*             &
                & 2/(Ppt%dx(5)+Ppt%dx(4))
        END IF

        IF (Ppt%StnclTplgy(6).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,6) = i
            Stencil%idx_clmns(MatStencil_j,6) = j+1
            Stencil%idx_clmns(MatStencil_k,6) = k
            Stencil%idx_val(6)=(Ppt%dx(4)*Ppt%dz(4))*                &
                & (Ppt%Cvt(6)%yy.ARMONIC.Ppt%Cvt(4)%yy)*             &
                & 2/(Ppt%dy(6)+Ppt%dy(4))
        END IF

        IF (Ppt%StnclTplgy(7).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,7) = i
            Stencil%idx_clmns(MatStencil_j,7) = j
            Stencil%idx_clmns(MatStencil_k,7) = k+1
            Stencil%idx_val(7)=(Ppt%dx(4)*Ppt%dy(4))*                &
                & (Ppt%Cvt(7)%zz.ARMONIC.Ppt%Cvt(4)%zz)*             &
                & 2/(Ppt%dz(7)+Ppt%dz(4))
        END IF

        Stencil%idx_val(4)= -(Stencil%idx_val(1)+Stencil%idx_val(2)+ &
            & Stencil%idx_val(3)+Stencil%idx_val(5)+                 &
            & Stencil%idx_val(6)+Stencil%idx_val(7))
        Stencil%idx_val(:)=Stencil%idx_val(:)/                       &
                         & (Ppt%dx(4)*Ppt%dy(4)*Ppt%dz(4))
        Stencil%IsActive=.TRUE.
        
    ELSEIF (Ppt%StnclTplgy(4).EQ.2) THEN ! Dirichlet cell
        Stencil%idx_val(4)=-one
    END IF


END SUBROUTINE GetTraditionalStencil

 !  - GetLiStencil: It's a routine that provides a Li stencil for a 
 !                  cell stored in a StensilVar.
 !    > IN: Ppt.
 !      + Ppt: It's a Property structure that contains every informa-
 !             tion needed to build the stencil on that cell.
 !    > OUT: Stencil, ierr.
 !      + Stencil: It's a StencilVar data structure that contains the 
 !                 information to be added to the matrix. See
 !                 StencilVar for more information.
 !      + ierr: It's an integer that indicates whether an error has 
 !              occurred during the call.

SUBROUTINE GetLiStencil(Ppt,Stencil,ierr)

    USE ANISOFLOW_Types,    ONLY : Property,StencilVar
    USE ANISOFLOW_Operators

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Property),INTENT(IN)               :: Ppt
    TYPE(StencilVar),INTENT(INOUT)          :: Stencil

    PetscInt                                :: i,j,k
    PetscReal                               :: one=1.d0,zero=0.d0

    IF (Ppt%StnclType.NE.2) THEN
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,               &
            & "ERROR: Properties must have a a stencil type box.\n"  &
            & ,ierr)
        STOP
    END IF

    ! The Li stencil is based on 19 cells.
    ALLOCATE(Stencil%idx_clmns(4,19))
    ALLOCATE(Stencil%idx_val(19))
    Stencil%idx_size=19

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
    Stencil%idx_val(:)=zero

    ! If the current cell is an active, source, or cauchy cell:
    IF ((Ppt%StnclTplgy(10).EQ.1).OR.(Ppt%StnclTplgy(10).EQ.3).OR.   &
        & (Ppt%StnclTplgy(10).EQ.4)) THEN

        ! 1-S Bloque izquierdo-centro-superior
        IF (Ppt%StnclTplgy(1).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,1) = i
            Stencil%idx_clmns(MatStencil_j,1) = j-1
            Stencil%idx_clmns(MatStencil_k,1) = k-1
            Stencil%idx_val(1)=                                      &
                & (Ppt%dy(10)*Ppt%dz(10))*                           &
                    & (Ppt%Cvt(7)%yz.ARMONIC.Ppt%Cvt(10)%yz)/        &
                    & (Ppt%dz(15)+2*Ppt%dz(7)+Ppt%dz(1))             &
                &+(Ppt%dx(10)*Ppt%dy(10))*                           &
                    & (Ppt%Cvt(3)%zy.ARMONIC.Ppt%Cvt(10)%zy)/        &
                    & (Ppt%dy(5 )+2*Ppt%dy(3)+Ppt%dy(1))
        END IF

        ! 2-O Bloque centro-detras-superior
        IF (Ppt%StnclTplgy(2).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,2) = i-1
            Stencil%idx_clmns(MatStencil_j,2) = j
            Stencil%idx_clmns(MatStencil_k,2) = k-1
            Stencil%idx_val(2)=                                      &
                & (Ppt%dy(10)*Ppt%dz(10))*                           &
                    & (Ppt%Cvt(9)%xz.ARMONIC.Ppt%Cvt(10)%xz)/        &
                    & (Ppt%dz(16)+2*Ppt%dz(9)+Ppt%dz(2))             &
                &+(Ppt%dx(10)*Ppt%dy(10))*                           &
                    & (Ppt%Cvt(3)%zx.ARMONIC.Ppt%Cvt(10)%zx)/        &
                    & (Ppt%dx(4 )+2*Ppt%dx(3)+Ppt%dx(2))
        END IF

        ! 3-J Bloque centro-centro-superior
        IF (Ppt%StnclTplgy(3).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,3) = i
            Stencil%idx_clmns(MatStencil_j,3) = j
            Stencil%idx_clmns(MatStencil_k,3) = k-1
            Stencil%idx_val(3)=                                      &
                & (Ppt%dy(10)*Ppt%dz(10))*                           &
                    & ((Ppt%Cvt(9)%xz.ARMONIC.Ppt%Cvt(10)%xz)-       &
                    & (Ppt%Cvt(11)%xz.ARMONIC.Ppt%Cvt(10)%xz))/      &
                    & (Ppt%dz(17)+2*Ppt%dz(10)+Ppt%dz(3))            &
                &+(Ppt%dx(10)*Ppt%dz(10))*                           &
                    & ((Ppt%Cvt(7)%yz.ARMONIC.Ppt%Cvt(10)%yz)-       &
                    & (Ppt%Cvt(13)%yz.ARMONIC.Ppt%Cvt(10)%yz))/      &
                    & (Ppt%dz(17)+2*Ppt%dz(10)+Ppt%dz(3))            &
                &+(Ppt%dx(10)*Ppt%dy(10))*                           &
                    & 2*(Ppt%Cvt(3)%zz.ARMONIC.Ppt%Cvt(10)%zz)/      &
                    & (Ppt%dz(3)+Ppt%dz(10))
        END IF

        ! 4-H Bloque centro-frontal-superior
        IF (Ppt%StnclTplgy(4).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,4) = i+1
            Stencil%idx_clmns(MatStencil_j,4) = j
            Stencil%idx_clmns(MatStencil_k,4) = k-1
            Stencil%idx_val(4)=                                      &
                & (Ppt%dy(10)*Ppt%dz(10))*                           &
                    & (-1)*(Ppt%Cvt(11)%xz.ARMONIC.Ppt%Cvt(10)%xz)/  &
                    & (Ppt%dz(18)+2*Ppt%dz(11)+Ppt%dz(4))            &
                &+(Ppt%dx(10)*Ppt%dy(10))*                           &
                    & (-1)*(Ppt%Cvt(3 )%zx.ARMONIC.Ppt%Cvt(10)%zx)/  &
                    & (Ppt%dx(4 )+2*Ppt%dx(3 )+Ppt%dx(2))
        END IF

        ! 5-Q Bloque derecho-centro-superior
        IF (Ppt%StnclTplgy(5).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,5) = i
            Stencil%idx_clmns(MatStencil_j,5) = j+1
            Stencil%idx_clmns(MatStencil_k,5) = k-1
            Stencil%idx_val(5)=                                      &
                & (Ppt%dy(10)*Ppt%dz(10))*                           &
                    & (-1)*(Ppt%Cvt(13)%yz.ARMONIC.Ppt%Cvt(10)%yz)/  &
                    & (Ppt%dz(19)+2*Ppt%dz(13)+Ppt%dz(5 ))           &
                &+(Ppt%dx(10)*Ppt%dy(10))*                           &
                    & (-1)*(Ppt%Cvt(3 )%zy.ARMONIC.Ppt%Cvt(10)%zy)/  &
                    & (Ppt%dy(5 )+2*Ppt%dy(3 )+Ppt%dy(1 ))
        END IF

        ! 6-M Bloque izquierdo-detras-centro
        IF (Ppt%StnclTplgy(6).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,6) = i-1
            Stencil%idx_clmns(MatStencil_j,6) = j-1
            Stencil%idx_clmns(MatStencil_k,6) = k
            Stencil%idx_val(6)=                                      &
                & (Ppt%dy(10)*Ppt%dz(10))*                           &
                    & (Ppt%Cvt(9)%xy.ARMONIC.Ppt%Cvt(10)%xy)/        &
                    & (Ppt%dy(12)+2*Ppt%dy(9)+Ppt%dy(6))             &
                &+(Ppt%dx(10)*Ppt%dz(10))*                           &
                    & (Ppt%Cvt(7)%yx.ARMONIC.Ppt%Cvt(10)%yx)/        &
                    & (Ppt%dx(8 )+2*Ppt%dx(6)+Ppt%dx(7))               ! subindices del diferencial malos en el paper(?) [6 <-> 7]
        END IF

        ! 7-F Bloque izquierdo-centro-centro
        IF (Ppt%StnclTplgy(7).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,7) = i
            Stencil%idx_clmns(MatStencil_j,7) = j-1
            Stencil%idx_clmns(MatStencil_k,7) = k
            Stencil%idx_val(7)=                                      &
                & (Ppt%dy(10)*Ppt%dz(10))*                           &
                    & ((Ppt%Cvt(9)%xy.ARMONIC.Ppt%Cvt(10)%xy)-       &
                    & (Ppt%Cvt(11)%xy.ARMONIC.Ppt%Cvt(10)%xy))/      &
                    & (Ppt%dy(13)+2*Ppt%dy(10)+Ppt%dy(7))            &
                &+(Ppt%dx(10)*Ppt%dz(10))*                           &
                    & 2*(Ppt%Cvt(7)%yy.ARMONIC.Ppt%Cvt(10)%yy)/      &
                    & (Ppt%dy(10)+Ppt%dy(7))                         &
                &+(Ppt%dx(10)*Ppt%dy(10))*                           &
                    & ((Ppt%Cvt(3)%zy.ARMONIC.Ppt%Cvt(10)%zy)-       &
                    & (Ppt%Cvt(17)%zy.ARMONIC.Ppt%Cvt(10)%zy))/      &
                    & (Ppt%dy(13)+2*Ppt%dy(10)+Ppt%dy(7))
        END IF

        ! 8-D Bloque izquierdo-frontal-centro
        IF (Ppt%StnclTplgy(8).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,8) = i+1
            Stencil%idx_clmns(MatStencil_j,8) = j-1
            Stencil%idx_clmns(MatStencil_k,8) = k
            Stencil%idx_val(8)=                                      &
                & (Ppt%dy(10)*Ppt%dz(10))*                           &
                    & (-1)*(Ppt%Cvt(11)%xy.ARMONIC.Ppt%Cvt(10)%xy)/  &
                    & (Ppt%dy(14)+2*Ppt%dy(8 )+Ppt%dy(11))           & ! subindices del diferencial malos en el paper(?) [8 <-> 11]
                &+(Ppt%dx(10)*Ppt%dz(10))*                           &
                    & (-1)*(Ppt%Cvt(7 )%yx.ARMONIC.Ppt%Cvt(10)%yx)/  &
                    & (Ppt%dx(8 )+2*Ppt%dx(6 )+Ppt%dx(7 ))             ! subindices del diferencial malos en el paper(?) [6 <-> 7]
        END IF

        ! 9-K Bloque centro-detras-centro
        IF (Ppt%StnclTplgy(9).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,9) = i-1
            Stencil%idx_clmns(MatStencil_j,9) = j
            Stencil%idx_clmns(MatStencil_k,9) = k
            Stencil%idx_val(9)=                                      &
                & (Ppt%dy(10)*Ppt%dz(10))*                           &
                    & 2*(Ppt%Cvt(9)%xx.ARMONIC.Ppt%Cvt(10)%xx)/      &
                    & (Ppt%dx(10)+Ppt%dx(9))                         &
                &+(Ppt%dx(10)*Ppt%dz(10))*                           &
                    & ((Ppt%Cvt(7)%yx.ARMONIC.Ppt%Cvt(10)%yx)-       &
                    & (Ppt%Cvt(13)%yx.ARMONIC.Ppt%Cvt(10)%yx))/      &
                    & (Ppt%dx(11)+2*Ppt%dx(10)+Ppt%dx(9))            &
                &+(Ppt%dx(10)*Ppt%dy(10))*                           &
                    & ((Ppt%Cvt(7)%zx.ARMONIC.Ppt%Cvt(10)%zx)-       &
                    & (Ppt%Cvt(13)%zx.ARMONIC.Ppt%Cvt(10)%zx))/      &
                    & (Ppt%dx(11)+2*Ppt%dx(10)+Ppt%dx(9))
        END IF

        ! 10-B Bloque centro-centro-centro
        Stencil%idx_val(10)=                                         &
            & (Ppt%dy(10)*Ppt%dz(10))*                               &
                & (-2)*( ((Ppt%Cvt(11)%xx.ARMONIC.Ppt%Cvt(10)%xx)/   &
                & (Ppt%dx(11)+Ppt%dx(10)))+                          &
                & ((Ppt%Cvt(9)%xx.ARMONIC.Ppt%Cvt(10)%xx)/           &
                & (Ppt%dx(10)+Ppt%dx(9))) )                          &
            &+(Ppt%dx(10)*Ppt%dz(10))*                               &
                & (-2)*( ((Ppt%Cvt(13)%yy.ARMONIC.Ppt%Cvt(10)%yy)/   &
                & (Ppt%dy(13)+Ppt%dy(10)))+                          &
                & ((Ppt%Cvt(7)%yy.ARMONIC.Ppt%Cvt(10)%yy)/           &
                & (Ppt%dy(10)+Ppt%dy(7))) )                          &
            &+(Ppt%dx(10)*Ppt%dy(10))*                               &
                & (-2)*( ((Ppt%Cvt(17)%zz.ARMONIC.Ppt%Cvt(10)%zz)/   &
                & (Ppt%dz(17)+Ppt%dz(10)))+                          &
                & ((Ppt%Cvt(3)%zz.ARMONIC.Ppt%Cvt(10)%zz)/           &
                & (Ppt%dz(10)+Ppt%dz(3))) )

        ! 11-A Bloque centro-frontal-centro
        IF (Ppt%StnclTplgy(11).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,11) = i+1
            Stencil%idx_clmns(MatStencil_j,11) = j
            Stencil%idx_clmns(MatStencil_k,11) = k
            Stencil%idx_val(11)=                                     &
                & (Ppt%dz(10)*Ppt%dy(10))*                           &
                    & 2*(Ppt%Cvt(11)%xx.ARMONIC.Ppt%Cvt(10)%xx)/     &
                    & (Ppt%dx(11)+Ppt%dx(10))                        &
                &+(Ppt%dx(10)*Ppt%dz(10))*                           &
                    & ((Ppt%Cvt(13)%yx.ARMONIC.Ppt%Cvt(10)%yx)-      &
                    & (Ppt%Cvt(7)%yx.ARMONIC.Ppt%Cvt(10)%yx))/       &
                    & (Ppt%dx(11)+2*Ppt%dx(10)+Ppt%dx(9))            &
                &+(Ppt%dx(10)*Ppt%dy(10))*                           &
                    & ((Ppt%Cvt(17)%zx.ARMONIC.Ppt%Cvt(10)%zx)-      &
                    & (Ppt%Cvt(3)%zx.ARMONIC.Ppt%Cvt(10)%zx))/       &
                    & (Ppt%dx(11)+2*Ppt%dx(10)+Ppt%dx(9))
        END IF

        ! 12-L Bloque derecho-detras-centro
        IF (Ppt%StnclTplgy(12).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,12) = i-1
            Stencil%idx_clmns(MatStencil_j,12) = j+1
            Stencil%idx_clmns(MatStencil_k,12) = k
            Stencil%idx_val(12)=                                     &
                & (Ppt%dy(10)*Ppt%dz(10))*                           &
                    & (-1)*(Ppt%Cvt(9 )%xy.ARMONIC.Ppt%Cvt(10)%xy)/  &
                    & (Ppt%dy(12)+2*Ppt%dy(9 )+Ppt%dy(6 ))           &
                &+(Ppt%dx(10)*Ppt%dz(10))*                           &
                    & (-1)*(Ppt%Cvt(13)%yx.ARMONIC.Ppt%Cvt(10)%yx)/  &
                    & (Ppt%dx(14)+2*Ppt%dx(12)+Ppt%dx(13))             ! subindices del diferencial malos en el paper(?) [12 <-> 13]
        END IF

        ! 13-E Bloque derecho-centro-centro
        IF (Ppt%StnclTplgy(13).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,13) = i
            Stencil%idx_clmns(MatStencil_j,13) = j+1
            Stencil%idx_clmns(MatStencil_k,13) = k
            Stencil%idx_val(13)=                                     &
                & (Ppt%dy(10)*Ppt%dz(10))*                           &
                    & ((Ppt%Cvt(11)%xy.ARMONIC.Ppt%Cvt(10)%xy)-      &
                    & (Ppt%Cvt(9 )%xy.ARMONIC.Ppt%Cvt(10)%xy))/      &
                    & (Ppt%dy(13)+2*Ppt%dy(10)+Ppt%dy(7))            &
                &+(Ppt%dx(10)*Ppt%dz(10))*                           &
                    & 2*(Ppt%Cvt(13)%yy.ARMONIC.Ppt%Cvt(10)%yy)/     &
                    & (Ppt%dy(13)+Ppt%dy(10))                        &
                &+(Ppt%dx(10)*Ppt%dy(10))*                           &
                    & ((Ppt%Cvt(17)%zy.ARMONIC.Ppt%Cvt(10)%zy)-      &
                    & (Ppt%Cvt(3 )%zy.ARMONIC.Ppt%Cvt(10)%zy))/      &
                    & (Ppt%dy(13)+2*Ppt%dy(10)+Ppt%dy(7))
        END IF

        ! 14-C Bloque derecho-frontal-centro
        IF (Ppt%StnclTplgy(14).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,14) = i+1
            Stencil%idx_clmns(MatStencil_j,14) = j+1
            Stencil%idx_clmns(MatStencil_k,14) = k
            Stencil%idx_val(14)=                                     &
                & (Ppt%dy(10)*Ppt%dz(10))*                           &
                    & (Ppt%Cvt(11)%xy.ARMONIC.Ppt%Cvt(10)%xy)/       &
                    & (Ppt%dy(14)+2*Ppt%dy(8 )+Ppt%dy(11))           & ! subindices del diferencial malos en el paper(?) [8  <-> 11]
                &+(Ppt%dx(10)*Ppt%dz(10))*                           &
                    & (Ppt%Cvt(13)%yx.ARMONIC.Ppt%Cvt(10)%yx)/       &
                    & (Ppt%dx(14)+2*Ppt%dx(12)+Ppt%dx(13))             ! subindices del diferencial malos en el paper(?) [12 <-> 13]
        END IF

        ! 15-R Bloque izquierdo-centro-inferior
        IF (Ppt%StnclTplgy(15).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,15) = i
            Stencil%idx_clmns(MatStencil_j,15) = j-1
            Stencil%idx_clmns(MatStencil_k,15) = k+1
            Stencil%idx_val(15)=                                     &
                & (Ppt%dy(10)*Ppt%dz(10))*                           &
                    & (-1)*(Ppt%Cvt(7 )%yz.ARMONIC.Ppt%Cvt(10)%yz)/  &
                    & (Ppt%dz(15)+2*Ppt%dz(7 )+Ppt%dz(1 ))           &
                &+(Ppt%dx(10)*Ppt%dy(10))*                           &
                    & (-1)*(Ppt%Cvt(17)%zy.ARMONIC.Ppt%Cvt(10)%zy)/  &
                    & (Ppt%dy(19)+2*Ppt%dy(17)+Ppt%dy(15))
        END IF

        ! 16-N Bloque centro-detras-inferior
        IF (Ppt%StnclTplgy(16).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,16) = i-1
            Stencil%idx_clmns(MatStencil_j,16) = j
            Stencil%idx_clmns(MatStencil_k,16) = k+1
            Stencil%idx_val(16)=                                     &
                & (Ppt%dy(10)*Ppt%dz(10))*                           &
                    & (-1)*(Ppt%Cvt(9 )%xz.ARMONIC.Ppt%Cvt(10)%xz)/  &
                    & (Ppt%dz(16)+2*Ppt%dz(9 )+Ppt%dz(2 ))           &
                &+(Ppt%dx(10)*Ppt%dy(10))*                           &
                    & (-1)*(Ppt%Cvt(17)%zx.ARMONIC.Ppt%Cvt(10)%zx)/  &
                    & (Ppt%dx(18)+2*Ppt%dx(17)+Ppt%dx(16))
        END IF

        ! 17-I Bloque centro-centro-inferior
        IF (Ppt%StnclTplgy(17).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,17) = i
            Stencil%idx_clmns(MatStencil_j,17) = j
            Stencil%idx_clmns(MatStencil_k,17) = k+1
            Stencil%idx_val(17)=                                     &
                & (Ppt%dy(10)*Ppt%dz(10))*                           &
                    & ((Ppt%Cvt(11)%xz.ARMONIC.Ppt%Cvt(10)%xz)-      &
                    & (Ppt%Cvt(9)%xz.ARMONIC.Ppt%Cvt(10)%xz))/       &
                    & (Ppt%dz(17)+2*Ppt%dz(10)+Ppt%dz(3))            &
                &+(Ppt%dx(10)*Ppt%dz(10))*                           &
                    & ((Ppt%Cvt(13)%yz.ARMONIC.Ppt%Cvt(10)%yz)-      &
                    & (Ppt%Cvt(7)%yz.ARMONIC.Ppt%Cvt(10)%yz))/       &
                    & (Ppt%dz(17)+2*Ppt%dz(10)+Ppt%dz(3))            &
                &+(Ppt%dx(10)*Ppt%dy(10))*                           &
                    & 2*(Ppt%Cvt(17)%zz.ARMONIC.Ppt%Cvt(10)%zz)/     &
                    & (Ppt%dz(17)+Ppt%dz(10))
        END IF

        ! 18-G Bloque centro-frontal-inferior
        IF (Ppt%StnclTplgy(18).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,18) = i+1
            Stencil%idx_clmns(MatStencil_j,18) = j
            Stencil%idx_clmns(MatStencil_k,18) = k+1
            Stencil%idx_val(18)=                                     &
                & (Ppt%dy(10)*Ppt%dz(10))*                           &
                    & (Ppt%Cvt(11)%xz.ARMONIC.Ppt%Cvt(10)%xz)/       &
                    & (Ppt%dz(18)+2*Ppt%dz(11)+Ppt%dz(4 ))           &
                &+(Ppt%dx(10)*Ppt%dy(10))*                           &
                    & (Ppt%Cvt(17)%zx.ARMONIC.Ppt%Cvt(10)%zx)/       &
                    & (Ppt%dx(18)+2*Ppt%dx(17)+Ppt%dx(16))
        END IF

        ! 19-P Bloque derecho-centro-inferior
        IF (Ppt%StnclTplgy(19).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,19) = i
            Stencil%idx_clmns(MatStencil_j,19) = j+1
            Stencil%idx_clmns(MatStencil_k,19) = k+1
            Stencil%idx_val(19)=                                     &
                & (Ppt%dy(10)*Ppt%dz(10))*                           &
                    & (Ppt%Cvt(13)%yz.ARMONIC.Ppt%Cvt(10)%yz)/       &
                    & (Ppt%dz(19)+2*Ppt%dz(13)+Ppt%dz(5 ))           &
                &+(Ppt%dx(10)*Ppt%dy(10))*                           &
                    & (Ppt%Cvt(17)%zy.ARMONIC.Ppt%Cvt(10)%zy)/       &
                    & (Ppt%dy(19)+2*Ppt%dy(17)+Ppt%dy(15))
        END IF

        Stencil%idx_val(:)=Stencil%idx_val(:)/                       &
                         & (Ppt%dx(10)*Ppt%dy(10)*Ppt%dz(10))
        Stencil%IsActive=.TRUE.

    ELSEIF (Ppt%StnclTplgy(10).EQ.2) THEN ! Dirichlet cell
        Stencil%idx_val(10)=-one
    END IF

END SUBROUTINE GetLiStencil


SUBROUTINE GetPerezStencil(Ppt,Stencil,ierr)

    USE ANISOFLOW_Types, ONLY : Property,StencilVar
    USE ANISOFLOW_Operators

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(Property),INTENT(IN)               :: Ppt
    TYPE(StencilVar),INTENT(INOUT)          :: Stencil

    PetscInt                                :: i,j,k
    PetscReal                               :: one=1.d0,zero=0.d0

    IF (Ppt%StnclType.NE.2) THEN
        CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,               &
            & "ERROR: Properties must have a a stencil type box.\n", &
            & ierr)
        STOP
    END IF

    ! The Li stencil is based on 19 cells.
    ALLOCATE(Stencil%idx_clmns(4,19))
    ALLOCATE(Stencil%idx_val(19))
    Stencil%idx_size=19

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
    Stencil%idx_val(:)=zero

    ! If the current cell is an active, source, or cauchy cell:
    IF ((Ppt%StnclTplgy(10).EQ.1).OR.(Ppt%StnclTplgy(10).EQ.3).OR.(Ppt%StnclTplgy(10).EQ.4)) THEN

        ! 1-P Bloque izquierdo-centro-superior
        IF (Ppt%StnclTplgy(1).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,1) = i
            Stencil%idx_clmns(MatStencil_j,1) = j-1
            Stencil%idx_clmns(MatStencil_k,1) = k-1
            Stencil%idx_val(1)=                                      &
                &-(Ppt%dx(10)*Ppt%dz(10))*                           &
                    & 2*(Ppt%Cvt(7)%yz.ARMONIC.Ppt%Cvt(10)%yz)/      &
                    & (Ppt%dz(15)+2*Ppt%dz(7)+Ppt%dz(1))             &
                &-(Ppt%dx(10)*Ppt%dy(10))*                           &
                    & 2*(Ppt%Cvt(3)%zy.ARMONIC.Ppt%Cvt(10)%zy)/      &
                    & (Ppt%dy(5 )+2*Ppt%dy(3)+Ppt%dy(1))              
        END IF

        ! 2-L Bloque centro-detras-superior
        IF (Ppt%StnclTplgy(2).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,2) = i-1
            Stencil%idx_clmns(MatStencil_j,2) = j
            Stencil%idx_clmns(MatStencil_k,2) = k-1
            Stencil%idx_val(2)=                                      &
                &-(Ppt%dy(10)*Ppt%dz(10))*                           &
                    & 2*(Ppt%Cvt(9)%xz.ARMONIC.Ppt%Cvt(10)%xz)/      &
                    & (Ppt%dz(16)+2*Ppt%dz(9)+Ppt%dz(2))             &
                &-(Ppt%dx(10)*Ppt%dy(10))*                           &
                    & 2*(Ppt%Cvt(3)%zx.ARMONIC.Ppt%Cvt(10)%zx)/      &
                    & (Ppt%dx(4 )+2*Ppt%dx(3)+Ppt%dx(2))
        END IF

        ! 3-F Bloque centro-centro-superior
        IF (Ppt%StnclTplgy(3).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,3) = i
            Stencil%idx_clmns(MatStencil_j,3) = j
            Stencil%idx_clmns(MatStencil_k,3) = k-1
            Stencil%idx_val(3)=                                      &
                & (Ppt%dx(10)*Ppt%dy(10))*                           &
                    & 2*(Ppt%Cvt(3)%zz.ARMONIC.Ppt%Cvt(10)%zz)/      &
                    & (Ppt%dz(3)+Ppt%dz(10))
        END IF

        ! 4-N Bloque centro-frontal-superior
        IF (Ppt%StnclTplgy(4).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,4) = i+1
            Stencil%idx_clmns(MatStencil_j,4) = j
            Stencil%idx_clmns(MatStencil_k,4) = k-1
            Stencil%idx_val(4)=                                      &
                & (Ppt%dy(10)*Ppt%dz(10))*                           &
                    & 2*(Ppt%Cvt(11)%xz.ARMONIC.Ppt%Cvt(10)%xz)/     &
                    & (Ppt%dz(18)+2*Ppt%dz(11)+Ppt%dz(4))            &
                &+(Ppt%dx(10)*Ppt%dy(10))*                           &
                    & 2*(Ppt%Cvt(3 )%zx.ARMONIC.Ppt%Cvt(10)%zx)/     &
                    & (Ppt%dx(4 )+2*Ppt%dx(3 )+Ppt%dx(2))
        END IF

        ! 5-R Bloque derecho-centro-superior
        IF (Ppt%StnclTplgy(5).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,5) = i
            Stencil%idx_clmns(MatStencil_j,5) = j+1
            Stencil%idx_clmns(MatStencil_k,5) = k-1
            Stencil%idx_val(5)=                                      &
                & (Ppt%dy(10)*Ppt%dz(10))*                           &
                    & 2*(Ppt%Cvt(13)%yz.ARMONIC.Ppt%Cvt(10)%yz)/     &
                    & (Ppt%dz(19)+2*Ppt%dz(13)+Ppt%dz(5 ))           &
                &+(Ppt%dx(10)*Ppt%dy(10))*                           &
                    & 2*(Ppt%Cvt(3 )%zy.ARMONIC.Ppt%Cvt(10)%zy)/     &
                    & (Ppt%dy(5 )+2*Ppt%dy(3 )+Ppt%dy(1 ))
        END IF

        ! 6-H Bloque izquierdo-detras-centro
        IF (Ppt%StnclTplgy(6).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,6) = i-1
            Stencil%idx_clmns(MatStencil_j,6) = j-1
            Stencil%idx_clmns(MatStencil_k,6) = k
            Stencil%idx_val(6)=                                      &
                &-(Ppt%dy(10)*Ppt%dz(10))*                           &
                    & 2*(Ppt%Cvt(9)%xy.ARMONIC.Ppt%Cvt(10)%xy)/      &
                    & (Ppt%dy(12)+2*Ppt%dy(9)+Ppt%dy(6))             &
                &-(Ppt%dx(10)*Ppt%dz(10))*                           &
                    & 2*(Ppt%Cvt(7)%xy.ARMONIC.Ppt%Cvt(10)%xy)/      &
                    & (Ppt%dx(6 )+2*Ppt%dx(7)+Ppt%dx(8))
        END IF

        ! 7-D Bloque izquierdo-centro-centro
        IF (Ppt%StnclTplgy(7).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,7) = i
            Stencil%idx_clmns(MatStencil_j,7) = j-1
            Stencil%idx_clmns(MatStencil_k,7) = k
            Stencil%idx_val(7)=                                      &
                & (Ppt%dx(10)*Ppt%dz(10))*                           &
                    & 2*(Ppt%Cvt(7)%yy.ARMONIC.Ppt%Cvt(10)%yy)/      &
                    & (Ppt%dy(10)+Ppt%dy(7))
        END IF

        ! 8-J Bloque izquierdo-frontal-centro
        IF (Ppt%StnclTplgy(8).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,8) = i+1
            Stencil%idx_clmns(MatStencil_j,8) = j-1
            Stencil%idx_clmns(MatStencil_k,8) = k
            Stencil%idx_val(8)=                                      &
                & (Ppt%dy(10)*Ppt%dz(10))*                           &
                    & 2*(Ppt%Cvt(11)%xy.ARMONIC.Ppt%Cvt(10)%xy)/     &
                    & (Ppt%dy(8)+2*Ppt%dy(11)+Ppt%dy(14))            &
                &+(Ppt%dx(10)*Ppt%dz(10))*                           &
                    & 2*(Ppt%Cvt(7 )%xy.ARMONIC.Ppt%Cvt(10)%xy)/     &
                    & (Ppt%dx(6)+2*Ppt%dx(7 )+Ppt%dx(8 ))
        END IF

        ! 9-B Bloque centro-detras-centro
        IF (Ppt%StnclTplgy(9).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,9) = i-1
            Stencil%idx_clmns(MatStencil_j,9) = j
            Stencil%idx_clmns(MatStencil_k,9) = k
            Stencil%idx_val(9)=                                      &
                & (Ppt%dy(10)*Ppt%dz(10))*                           &
                    & 2*(Ppt%Cvt(9)%xx.ARMONIC.Ppt%Cvt(10)%xx)/      &
                    & (Ppt%dx(10)+Ppt%dx(9))
        END IF

        ! 11-C Bloque centro-frontal-centro
        IF (Ppt%StnclTplgy(11).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,11) = i+1
            Stencil%idx_clmns(MatStencil_j,11) = j
            Stencil%idx_clmns(MatStencil_k,11) = k
            Stencil%idx_val(11)=                                     &
                & (Ppt%dz(10)*Ppt%dy(10))*                           &
                    & 2*(Ppt%Cvt(11)%xx.ARMONIC.Ppt%Cvt(10)%xx)/     &
                    & (Ppt%dx(11)+Ppt%dx(10))
        END IF

        ! 12-I Bloque derecho-detras-centro
        IF (Ppt%StnclTplgy(12).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,12) = i-1
            Stencil%idx_clmns(MatStencil_j,12) = j+1
            Stencil%idx_clmns(MatStencil_k,12) = k
            Stencil%idx_val(12)=                                     &
                & (Ppt%dy(10)*Ppt%dz(10))*                           &
                    & 2*(Ppt%Cvt(9 )%xy.ARMONIC.Ppt%Cvt(10)%xy)/     &
                    & (Ppt%dy(6 )+2*Ppt%dy(9 )+Ppt%dy(12))           &
            &+(Ppt%dx(10)*Ppt%dz(10))*                               &
                & 2*(Ppt%Cvt(13)%xy.ARMONIC.Ppt%Cvt(10)%xy)/         &
                & (Ppt%dx(12)+2*Ppt%dx(13)+Ppt%dx(14))
        END IF

        ! 13-E Bloque derecho-centro-centro
        IF (Ppt%StnclTplgy(13).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,13) = i
            Stencil%idx_clmns(MatStencil_j,13) = j+1
            Stencil%idx_clmns(MatStencil_k,13) = k
            Stencil%idx_val(13)=                                     &
                & (Ppt%dx(10)*Ppt%dz(10))*                           &
                    & 2*(Ppt%Cvt(13)%yy.ARMONIC.Ppt%Cvt(10)%yy)/     &
                    & (Ppt%dy(13)+Ppt%dy(10))
        END IF

        ! 14-K Bloque derecho-frontal-centro
        IF (Ppt%StnclTplgy(14).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,14) = i+1
            Stencil%idx_clmns(MatStencil_j,14) = j+1
            Stencil%idx_clmns(MatStencil_k,14) = k
            Stencil%idx_val(14)=                                     &
                &-(Ppt%dy(10)*Ppt%dz(10))*                           &
                    & 2*(Ppt%Cvt(11)%xy.ARMONIC.Ppt%Cvt(10)%xy)/     &
                    & (Ppt%dy(8 )+2*Ppt%dy(11)+Ppt%dy(14))           &
                &-(Ppt%dx(10)*Ppt%dz(10))*                           &
                    & 2*(Ppt%Cvt(13)%xy.ARMONIC.Ppt%Cvt(10)%xy)/     &
                    & (Ppt%dx(12)+2*Ppt%dx(13)+Ppt%dx(14))
        END IF

        ! 15-Q Bloque izquierdo-centro-inferior
        IF (Ppt%StnclTplgy(15).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,15) = i
            Stencil%idx_clmns(MatStencil_j,15) = j-1
            Stencil%idx_clmns(MatStencil_k,15) = k+1
            Stencil%idx_val(15)=                                     &
                & (Ppt%dx(10)*Ppt%dz(10))*                           &
                    & 2*(Ppt%Cvt(7 )%yz.ARMONIC.Ppt%Cvt(10)%yz)/     &
                    & (Ppt%dz(15)+2*Ppt%dz(7 )+Ppt%dz(1 ))           &
                &+(Ppt%dx(10)*Ppt%dy(10))*                           &
                    & 2*(Ppt%Cvt(17)%zy.ARMONIC.Ppt%Cvt(10)%zy)/     &
                    & (Ppt%dy(19)+2*Ppt%dy(17)+Ppt%dy(15))
        END IF

        ! 16-M Bloque centro-detras-inferior
        IF (Ppt%StnclTplgy(16).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,16) = i-1
            Stencil%idx_clmns(MatStencil_j,16) = j
            Stencil%idx_clmns(MatStencil_k,16) = k+1
            Stencil%idx_val(16)=                                     &
                & (Ppt%dy(10)*Ppt%dz(10))*                           &
                    & 2*(Ppt%Cvt(9 )%xz.ARMONIC.Ppt%Cvt(10)%xz)/     &
                    & (Ppt%dz(16)+2*Ppt%dz(9 )+Ppt%dz(2 ))           &
                &+(Ppt%dx(10)*Ppt%dy(10))*                           &
                    & 2*(Ppt%Cvt(17)%zx.ARMONIC.Ppt%Cvt(10)%zx)/     &
                    & (Ppt%dx(18)+2*Ppt%dx(17)+Ppt%dx(16))
        END IF

        ! 17-I Bloque centro-centro-inferior
        IF (Ppt%StnclTplgy(17).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,17) = i
            Stencil%idx_clmns(MatStencil_j,17) = j
            Stencil%idx_clmns(MatStencil_k,17) = k+1
            Stencil%idx_val(17)=                                     &
                & (Ppt%dx(10)*Ppt%dy(10))*                           &
                    & 2*(Ppt%Cvt(17)%zz.ARMONIC.Ppt%Cvt(10)%zz)/     &
                    & (Ppt%dz(17)+Ppt%dz(10))
        END IF

        ! 18-O Bloque centro-frontal-inferior
        IF (Ppt%StnclTplgy(18).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,18) = i+1
            Stencil%idx_clmns(MatStencil_j,18) = j
            Stencil%idx_clmns(MatStencil_k,18) = k+1
            Stencil%idx_val(18)=                                     &
                &-(Ppt%dy(10)*Ppt%dz(10))*                           &
                    & 2*(Ppt%Cvt(11)%xz.ARMONIC.Ppt%Cvt(10)%xz)/     &
                    & (Ppt%dz(18)+2*Ppt%dz(11)+Ppt%dz(4 ))           &
                &-(Ppt%dx(10)*Ppt%dy(10))*                           &
                    & 2*(Ppt%Cvt(17)%zx.ARMONIC.Ppt%Cvt(10)%zx)/     &
                    & (Ppt%dx(18)+2*Ppt%dx(17)+Ppt%dx(16))
        END IF

        ! 19-S Bloque derecho-centro-inferior
        IF (Ppt%StnclTplgy(19).NE.0) THEN
            Stencil%idx_clmns(MatStencil_i,19) = i
            Stencil%idx_clmns(MatStencil_j,19) = j+1
            Stencil%idx_clmns(MatStencil_k,19) = k+1
            Stencil%idx_val(19)=                                     &
                &-(Ppt%dx(10)*Ppt%dz(10))*                           &
                    & 2*(Ppt%Cvt(13)%yz.ARMONIC.Ppt%Cvt(10)%yz)/     &
                    & (Ppt%dz(15)+2*Ppt%dz(13)+Ppt%dz(5 ))           &
                &-(Ppt%dx(10)*Ppt%dy(10))*                           &
                    & 2*(Ppt%Cvt(3 )%zy.ARMONIC.Ppt%Cvt(10)%zy)/     &
                    & (Ppt%dy(19)+2*Ppt%dy(17)+Ppt%dy(15))              
        END IF

        ! 10-A Bloque centro-centro-centro
        Stencil%idx_val(10)= -(Stencil%idx_val(3)+Stencil%idx_val(7)+&
            & Stencil%idx_val(9)+Stencil%idx_val(11)+                &
            & Stencil%idx_val(13)+Stencil%idx_val(17))
        Stencil%idx_val(:)=Stencil%idx_val(:)/                       &
                         & (Ppt%dx(10)*Ppt%dy(10)*Ppt%dz(10))
        Stencil%IsActive=.TRUE.

    ELSEIF (Ppt%StnclTplgy(10).EQ.2) THEN ! Dirichlet cell
        Stencil%idx_val(10)=-one
    END IF

END SUBROUTINE GetPerezStencil

SUBROUTINE ApplyDirichlet(BCFld,TimeZone,A,b,ierr)

    USE ANISOFLOW_Types,     ONLY : Geometry,BoundaryConditions
    USE ANISOFLOW_Interface, ONLY : GetVerbose

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(BoundaryConditions),INTENT(IN)     :: BCFld
    PetscInt,INTENT(IN)                     :: TimeZone
    Mat,INTENT(INOUT)                       :: A
    Vec,INTENT(INOUT)                       :: b

    VecScatter                              :: Scatter
    IS                                      :: NaturalOrder
    PetscInt                                :: DirichSize
    PetscReal                               :: one=1.D0
    PetscBool                               :: Verbose


    CALL MatZeroRowsIS(A,BCFld%DirichIS(TimeZone),-one,              &
        & PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)

    CALL VecGetSize(BCFld%Dirich(TimeZone),DirichSize,ierr)

    CALL ISCreateStride(PETSC_COMM_WORLD,DirichSize,0,1,NaturalOrder,& 
        & ierr)
    CALL VecScatterCreate(BCFld%Dirich(TimeZone),NaturalOrder,b,     & 
        & BCFld%DirichIS(TimeZone),Scatter,ierr)

    CALL VecScatterBegin(Scatter,BCFld%Dirich(TimeZone),b,           &
        & INSERT_VALUES,SCATTER_FORWARD,ierr)
    CALL VecScatterEnd(Scatter,BCFld%Dirich(TimeZone),b,             &
        & INSERT_VALUES,SCATTER_FORWARD,ierr)

    CALL VecScatterDestroy(Scatter,ierr)
    CALL ISDestroy(NaturalOrder,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,      &
        & "[GetSystem Event] Dirichlet boundary conditions properl"//&
        & "y implemented\n",ierr)

END SUBROUTINE ApplyDirichlet

SUBROUTINE ApplySource(BCFld,TimeZone,b,ierr)

    USE ANISOFLOW_Types,     ONLY : Geometry,BoundaryConditions
    USE ANISOFLOW_Interface, ONLY : GetVerbose

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

    CALL ISCreateStride(PETSC_COMM_WORLD,SourceSize,0,1,NaturalOrder,&
        & ierr)
    CALL VecScatterCreate(BCFld%Source(TimeZone),NaturalOrder,b,     &
        & BCFld%SourceIS(TimeZone),Scatter,ierr)

    CALL VecScatterBegin(Scatter,BCFld%Source(TimeZone),b,ADD_VALUES,&
        & SCATTER_FORWARD,ierr)
    CALL VecScatterEnd(Scatter,BCFld%Source(TimeZone),b,ADD_VALUES,  &
        & SCATTER_FORWARD,ierr)

    CALL VecScatterDestroy(Scatter,ierr)
    CALL ISDestroy(NaturalOrder,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,      &
        & "[GetSystem Event] Source terms properly implemented\n",   &
        & ierr)

END SUBROUTINE ApplySource

SUBROUTINE ApplyCauchy(BCFld,TimeZone,A,b,ierr)

    USE ANISOFLOW_Types,     ONLY : Geometry,BoundaryConditions
    USE ANISOFLOW_Interface, ONLY : GetVerbose

    IMPLICIT NONE

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    PetscErrorCode,INTENT(INOUT)            :: ierr
    TYPE(BoundaryConditions),INTENT(IN)     :: BCFld
    PetscInt,INTENT(IN)                     :: TimeZone
    Vec,INTENT(INOUT)                       :: b
    Mat,INTENT(INOUT)                       :: A

    Vec                                     :: MultHC,MultHeC,DiagA, &
                                             & DiagApartial
    VecScatter                              :: Scatter
    IS                                      :: NaturalOrder
    PetscBool                               :: Verbose

!     CALL VecView(BCFld%CauchyHe(TimeZone),PETSC_VIEWER_STDOUT_WORLD, &
!         & ierr)
!     CALL VecView(BCFld%CauchyC(TimeZone),PETSC_VIEWER_STDOUT_WORLD,  &
!         & ierr)

    CALL VecDuplicate(BCFld%CauchyHe(TimeZone),MultHC,ierr)
    CALL VecDuplicate(BCFld%CauchyHe(TimeZone),MultHeC,ierr)
    CALL VecDuplicate(BCFld%CauchyHe(TimeZone),DiagApartial,ierr)

    CALL VecPointwiseMult(MultHeC,BCFld%CauchyC(TimeZone),           &
        & BCFld%CauchyHe(TimeZone),ierr)

!     CALL VecView(MultHeC,PETSC_VIEWER_STDOUT_WORLD,ierr)

    CALL ISCreateStride(PETSC_COMM_WORLD,BCFld%SizeCauchy,0,1,       &
        & NaturalOrder,ierr)
    CALL VecScatterCreate(MultHeC,NaturalOrder,b,                    &
        & BCFld%CauchyIS(TimeZone),Scatter,ierr)

    CALL VecScatterBegin(Scatter,MultHeC,b,ADD_VALUES,               &
        & SCATTER_FORWARD,ierr)
    CALL VecScatterEnd(Scatter,MultHeC,b,ADD_VALUES,SCATTER_FORWARD, &
        & ierr)

    CALL VecDestroy(MultHeC,ierr)

    CALL VecDuplicate(b,DiagA,ierr)
    CALL MatGetDiagonal(A,DiagA,ierr)

    CALL VecScatterBegin(Scatter,DiagA,DiagApartial,ADD_VALUES,      &
        & SCATTER_REVERSE,ierr)
    CALL VecScatterEnd(Scatter,DiagA,DiagApartial,ADD_VALUES,        &
        & SCATTER_REVERSE,ierr)

!     CALL VecView(DiagA,PETSC_VIEWER_STDOUT_WORLD,ierr)
!     CALL VecView(DiagApartial,PETSC_VIEWER_STDOUT_WORLD,ierr)

    CALL VecPointwiseMult(DiagApartial,BCFld%CauchyC(TimeZone),MultHC,&
        & ierr)

    CALL VecScatterBegin(Scatter,DiagApartial,DiagA,ADD_VALUES,      &
        & SCATTER_FORWARD,ierr)
    CALL VecScatterEnd(Scatter,DiagApartial,DiagA,ADD_VALUES,        &
        & SCATTER_FORWARD,ierr)

    CALL MatDiagonalSet(A,diagA,INSERT_VALUES,ierr)

    CALL VecScatterDestroy(Scatter,ierr)
    CALL ISDestroy(NaturalOrder,ierr)

    CALL GetVerbose(Verbose,ierr)
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,      &
        & "[GetSystem Event] Cauchy terms properly implemented\n",ierr)

END SUBROUTINE ApplyCauchy

SUBROUTINE ApplyTimeDiff(PptFld,BCFld,TimeZone,TimeStep,A,b,x,ierr)

    USE ANISOFLOW_Types,     ONLY : PropertiesField,BoundaryConditions
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

    Vec                                     :: VecDT,VecZero
    PetscReal                               :: DT,one=1.d0,zero=0.d0
    CHARACTER(LEN=200)                      :: CharDT
    PetscBool                               :: Verbose
    PetscInt                                :: DirichSize
    IS                                      :: NaturalOrder
    VecScatter                              :: Scatter


    CALL GetVerbose(Verbose,ierr)
    CALL GetDT(BCFld,TimeZone,TimeStep,DT,ierr)
    WRITE(CharDT,*)DT
    IF (Verbose) CALL PetscSynchronizedPrintf(PETSC_COMM_WORLD,      &
        & "[GetSystem Stage] DT: "//ADJUSTL(TRIM(CharDT))//"\n",ierr)

    CALL VecDuplicate(x,VecDT,ierr)
    CALL VecSet(VecDT,-one/DT,ierr)
    CALL VecPointwiseMult(VecDT,VecDT,PptFld%Sto%Cell,ierr)

    CALL VecGetSize(BCFld%Dirich(TimeZone),DirichSize,ierr)

    ! Inserting zeros to maintain dirichlet BC
    CALL VecDuplicate(BCFld%Dirich(TimeZone),VecZero,ierr)
    CALL VecSet(VecZero,zero,ierr)
    CALL ISCreateStride(PETSC_COMM_WORLD,DirichSize,0,1,NaturalOrder,&
        & ierr)
    CALL VecScatterCreate(VecZero,NaturalOrder,VecDT,                &
        & BCFld%DirichIS(TimeZone),Scatter,ierr)
    CALL VecScatterBegin(Scatter,VecZero,VecDT,INSERT_VALUES,        &
        & SCATTER_FORWARD,ierr)
    CALL VecScatterEnd(Scatter,VecZero,VecDT,INSERT_VALUES,          &
        & SCATTER_FORWARD,ierr)

    CALL MatDiagonalSet(A,VecDT,ADD_VALUES,ierr)
    CALL VecPointwiseMult(VecDT,VecDT,x,ierr)
    CALL VecAXPY(b,one,VecDT,ierr)

    CALL VecScatterDestroy(Scatter,ierr)
    CALL ISDestroy(NaturalOrder,ierr)
    CALL VecDestroy(VecDT,ierr)
    CALL VecDestroy(VecZero,ierr)


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