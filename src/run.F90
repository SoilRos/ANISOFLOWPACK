!ORDEN DE EJECUCION:
!Creacion de variable tipo geometria
!Preproceasmiento de la informacion segun sea para reimgen transitorio o permanente
!Creacion de topologia
PROGRAM run
    !
    USE anisotype
    USE geometry
    USE iodata
    USE topology
    USE isosystem
    USE anisosystem
    USE k_operations
    !
    IMPLICIT NONE	
#include <petsc/finclude/petscsys.h>
	!Declaracion de variables de petsc
	PetscErrorCode								:: ierr
	PetscMPIInt									:: tprocess, process
    !Declaracion de variables tipo
    TYPE(domain)                                :: dmngmtry
    TYPE(topology)                              :: topo
    !Declaracion de variables tipo caracter
    CHARACTER(200)                              :: prjctrt,prjctrtd,route,rt
    CHARACTER(20)                               :: flnm,forminI,forminR
    CHARACTER(1) :: intro=' '  
    CHARACTER(8)::prjctnm
    !Declaracion de constantes
    INTEGER(I4B)                                :: i,j,k,cont,vl,aux,u,n
    REAL(DP)                                    :: xi,yi,zi,xf,yf,zf,dx,dy,dz
    !Declaracion de arreglos
    REAL(DP), DIMENSION(:), ALLOCATABLE         :: X,Y,Z,kxx,kyy,kzz,kxy,kxz,kyz,qhext
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE     :: topheconect
    INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE   :: ijheconect
    REAL(DP), DIMENSION(:,:), ALLOCATABLE       :: Teq,qhcte
    !    
    PARAMETER(u=10)
    !Anisotropia 10
    PARAMETER(prjctrt='/cygdrive/D/Colciencias-Tesis/Validacion/Ans_10/Project')
    PARAMETER(prjctrtd='/cygdrive/D/Colciencias-Tesis/Validacion/Ans_10/Data')
!    !Anisotropia 20
!    PARAMETER(prjctrt='/cygdrive/D/Colciencias-Tesis/Validacion/Ans_40/Project')
!    PARAMETER(prjctrtd='/cygdrive/D/Colciencias-Tesis/Validacion/Ans_40/Data')
    
!    !Anisotropia 30
!    PARAMETER(prjctrt='/cygdrive/D/Colciencias-Tesis/Validacion/Ans_40/Project')
!    PARAMETER(prjctrtd='/cygdrive/D/Colciencias-Tesis/Validacion/Ans_40/Data')
    
!    !Anisotropia 40
!    PARAMETER(prjctrt='/cygdrive/D/Colciencias-Tesis/Validacion/Ans_40/Project')
!    PARAMETER(prjctrtd='/cygdrive/D/Colciencias-Tesis/Validacion/Ans_40/Data')
    
!    !isotropia
!    PARAMETER(prjctrt='/cygdrive/D/Colciencias-Tesis/Validacion/Isotropy/Project')
!    PARAMETER(prjctrtd='/cygdrive/D/Colciencias-Tesis/Validacion/Isotropy/Data')


	call PetscInitialize(PETSC_COMM_WORLD,ierr) 
	call MPI_comm_rank(MPI_COMM_WORLD,process,ierr)
	call MPI_comm_size(MPI_COMM_WORLD,tprocess,ierr)

	if (process==0) then


	

    flnm='sanpck.domnRST' 
    forminI='I6'
    forminR='ES14.8'

    CALL bld_domainff(prjctrtd,flnm,X,Y,Z,dmngmtry) 
!    CALL imprtdtstbn(prjctrtd,prjctrt,forminI,forminR,dmngmtry)
    CALL aniso_imprtdtstbn(prjctrtd,prjctrt,forminI,forminR,dmngmtry)
    CALL K_tensor(prjctrt,dmngmtry,kxx,kyy,kzz,kxy,kxz,kyz)
!    CALL bld_topologybn(prjctrt,dmngmtry,topo)   
!    CALL isoTdns_bn(prjctrt,topo,dmngmtry,Teq,topheconect,ijheconect,qhcte,qhext)
!    aux=dmngmtry%clmns*dmngmtry%rws*dmngmtry%lvls
!    rt=trim(prjctrt)//'/sanpck.alf'
!    OPEN(u, FILE=rt, STATUS='OLD', ACCESS='DIRECT', RECL=18*8, FORM='UNFORMATTED')
!    CALL gtvlr_bn(u,dmngmtry,2,2,2,vl) 
    
!   Bucle para imprimir matriz densa
!    do i=topo%act/2,topo%act
!        do j=topo%act/2,topo%act
!            write(*,*)Teq(i,j),i,j
!        end do
!    end do
    
!    UTIL PARA ABRIR UN ARCHIVO Y LEER
!    rt=trim(prjctrt)//'/sanpck.kxx'
!    OPEN(u, FILE=rt, STATUS='OLD', ACCESS='DIRECT', RECL=18*8, FORM='UNFORMATTED')
!    n=dmngmtry%lvls*dmngmtry%rws*dmngmtry%clmns
!    DO i=1,n
!        READ(u,rec=i)vl
!        WRITE(*,*)vl
!    END DO

	end if
CALL Petscfinalize(ierr)

END PROGRAM

