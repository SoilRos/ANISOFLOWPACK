!-----------------------------------------------------------------------------------------------------------------------------------
!Copyright (c) 2013, Alvarez-Villa O.D., Perez Tevin  GOTTA Ingenieria SAS, All Right Reserved.
!
!This Program is distributed in the hope that they will be used but WITHOUT ANY WARRANTY. No autor or distributor accepts 
!responsability to anyone for consequences of using them or for whether they serve !any particular porpouse or work at all, unless he
!says so in writing. Everyone is granted permission to copy, modify an distribute this program, but only under the conditions that 
!this notice and above copyright notice remains intact.
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!Este modulo contiene las variables tipo usadas en ANISOFLOWPACK. Las variables tipo implementadas han sido adaptadas partir de los
!modulos de variables tipo de los siguientes trabajos:
!1. Festellustype del software de interpolacion espacial usando tecnicas geoestadisticas Festellus. 2014, Rendon-Alvarez J.P.
!2. Eigentype del software para modelar flujo subterraneo FDPACK. 2013, Alvarez-Villa O.D.
!-----------------------------------------------------------------------------------------------------------------------------------
!
MODULE isosystem
    INTERFACE ensmblT_bn
        SUBROUTINE isoTdns_bn(prjctrt,topo,dmngmtry,Teq,topheconect,ijheconect,qhcte,qhext)
            USE anisotype
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN)                                       :: prjctrt
            TYPE(topology), INTENT(IN)                                     :: topo
            TYPE(domain), INTENT(IN)                                       :: dmngmtry
            REAL(DP), DIMENSION(:,:),ALLOCATABLE, INTENT(OUT)              :: Teq
            INTEGER(I4B), DIMENSION(:),ALLOCATABLE, INTENT(OUT)            :: topheconect
            INTEGER(I4B), DIMENSION(:,:),ALLOCATABLE, INTENT(OUT)          :: ijheconect
            REAL(DP), DIMENSION(:,:),ALLOCATABLE, INTENT(OUT)              :: qhcte
            REAL(DP), DIMENSION(:),ALLOCATABLE, INTENT(OUT)                :: qhext
        END SUBROUTINE isoTdns_bn
        !BL
        SUBROUTINE isoTdns_bn_petsc(prjctrt,topo,dmngmtry,Teq,topheconect,ijheconect,qhcte,qhext)
            USE anisotype
            IMPLICIT NONE
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscvec.h>
            CHARACTER(*), INTENT(IN)                                       :: prjctrt
            TYPE(topology), INTENT(IN)                                     :: topo
            TYPE(domain), INTENT(IN)                                       :: dmngmtry
            Mat, INTENT(OUT)                                               :: Teq
            INTEGER(I4B), DIMENSION(:),ALLOCATABLE, INTENT(OUT)            :: topheconect
            INTEGER(I4B), DIMENSION(:,:),ALLOCATABLE, INTENT(OUT)          :: ijheconect
            REAL(DP), DIMENSION(:,:),ALLOCATABLE, INTENT(OUT)              :: qhcte
            REAL(DP), DIMENSION(:),ALLOCATABLE, INTENT(OUT)                :: qhext
        END SUBROUTINE isoTdns_bn_petsc

        SUBROUTINE isoTsprc_bn(prjctrt,topo,dmngmtry,steq,topheconect,ijheconect,qhcte,qhext)
            USE anisotype
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN)                           :: prjctrt
            TYPE(topology), INTENT(IN)                         :: topo
            TYPE(domain), INTENT(IN)                           :: dmngmtry
            TYPE(dispersed), INTENT(OUT)                       :: steq
            INTEGER(I4B), DIMENSION(:), INTENT(OUT)            :: topheconect
            INTEGER(I4B), DIMENSION(:,:), INTENT(OUT)          :: ijheconect
            REAL(DP), DIMENSION(:,:), INTENT(OUT)              :: qhcte
            REAL(DP), DIMENSION(:), INTENT(OUT)                :: qhext
        END SUBROUTINE isoTsprc_bn
    END INTERFACE ensmblT_bn
    !!
    INTERFACE
        SUBROUTINE isoSF(prjctrt,topo,dmngmtry,SF)
            USE anisotype
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN)                           :: prjctrt
            TYPE(topology), INTENT(IN)                         :: topo
            TYPE(domain), INTENT(IN)                           :: dmngmtry
            REAL(DP), DIMENSION(:), INTENT(INOUT)              :: SF
        END SUBROUTINE isoSF
    END INTERFACE
    !!
    INTERFACE
        SUBROUTINE isochangehcterow(p,stts,qub,nheconect,topheconect,ijheconect,qhcte,a,i,j,pstn)
            USE anisotype
            IMPLICIT NONE
            INTEGER(I4B), INTENT(INOUT)                          :: p
            INTEGER(I4B), INTENT(IN)                             :: stts
            REAL(DP), INTENT(IN)                                 :: qub
            INTEGER(I4B), INTENT(INOUT)                          :: nheconect
            INTEGER(I4B), DIMENSION(:), INTENT(INOUT)            :: topheconect
            INTEGER(I4B), DIMENSION(:,:), INTENT(INOUT)          :: ijheconect
            REAL(DP), DIMENSION(:,:), INTENT(INOUT)              :: qhcte
            INTEGER(I4B), INTENT(IN)                             :: a,i,j,pstn
        END SUBROUTINE isochangehcterow
    END INTERFACE
!    !!
END MODULE isosystem
!
!***********************************************************************************************************************************
!*******************SUBRUTINA QUE ENSAMBLA LA MATRIZ DE TRANSMISIVIDADES EQUIVALENTES EN NUMEROS BINARIOS***************************
!***************************@author Álvarez-Villa O.D. & Perez K.*******************************************************************
!*********************************************************************************************************************************** 
!
!DESCRIPCION:
!NOMBRE: isoTdns_bn

!Subrutina para ensamblar la matriz de transmisividades equivalentes a partir de las tansmisividades de bloque y
!el vector de flujos unitarios generados por los bloques de altura piezomotrica prescita sobre los bloques activos.
!
!VARIABLES GLOBALES:
!ENTRADA:
!prjctrt        :Ruta del proyecto de trabajo.
!topo           :La topologia del acuifero
!dmngmtry       :La geometria del acuifero
!SALIDAS:
!Teq            :La matriz de transmisividades equivalentes
!topheconect    :Arreglo de indicadores de celdas activas donde existe conexion.
!ijhconect      :Arreglo que indican cuales celdas con altura impuesta se conectan con las activas.
!qhcte          :Arreglo con las conductividades equivalentes de las celdas activas con conexion.
!qhext          :Flujos generados por las alturas externas sobre las celdas activas
!
!VARIBALES LOCALES:
!Ka= conductividad entre bloque central e i-1
!Kb= conductividad entre bloque central e i+1
!Kc= conductividad entre bloque central e j-1
!Kd= conductividad entre bloque central e j-1
!Ke= conductividad entre bloque central e j+1
!Kf= conductividad entre bloque central e k-1
!Kg= conductividad entre bloque central e k+1
!Kr= conductividad entre bloque central y celda de nivel externo.
!E=Espesor de celda de nivel externo
!Le=Longitud de celda de nivel externo
!Wd=Ancho de celda de nivel externo
!dx,dy,dz=longitud del bloque central en cada direccion cartesiana
!qub=Variable auxiliar para acumular flujo
!dist=Vairiable auxiliar que almacena distancia entre bloques
!armncmn=Varibale auxiliar que almacena media armonica
!

SUBROUTINE isoTdns_bn(prjctrt,topo,dmngmtry,Teq,topheconect,ijheconect,qhcte,qhext)
    USE anisotype
    USE geometry
    USE iodata,     ONLY : gtvlr_bn
    USE isosystem,  ONLY : isochangehcterow
    USE util
    !
    IMPLICIT NONE
    !
    CHARACTER(*), INTENT(IN)                                       :: prjctrt
    TYPE(topology), INTENT(IN)                                     :: topo
    TYPE(domain), INTENT(IN)                                       :: dmngmtry
    REAL(DP), DIMENSION(:,:),ALLOCATABLE, INTENT(OUT)              :: Teq
    INTEGER(I4B), DIMENSION(:),ALLOCATABLE, INTENT(OUT)            :: topheconect
    INTEGER(I4B), DIMENSION(:,:),ALLOCATABLE, INTENT(OUT)          :: ijheconect
    REAL(DP), DIMENSION(:,:),ALLOCATABLE, INTENT(OUT)              :: qhcte
    REAL(DP), DIMENSION(:),ALLOCATABLE, INTENT(OUT)                :: qhext
    !
    REAL(DP)       :: Ka,Kb,Kc,Kd,Ke,Kf,Kg,Kr,E,Le,Wd,dx,dy,dz,qub,dist
    CHARACTER(200) :: rtact,rthe,rttrns,rtcndtnc,rtdpth,rtl,rtb
    INTEGER(I4B)   :: a,b,i,j,k,m,l,o,u,v,w,ua,ul,ub,uhe,stts,sthe,n,nheconect,p,dim
    !
    PARAMETER(ua=37);  PARAMETER(u=65)
    PARAMETER(v=58);   PARAMETER(w=33)
    PARAMETER(ul=38);  PARAMETER(ub=39)
    PARAMETER(uhe=74)
    nheconect=0; p=0; 
    !
    WRITE (*,*) '*---|Inicia construccion matriz densa de conductividades|---*'
    !
    !Revisa si la topologia y geometria se han construido
    IF(dmngmtry%exist.EQV..false.)CALL error('En ensmblT: La geometria no se ha construido')
    IF(topo%hastopo.EQV..false.)CALL error('En ensmblT: La topologia no se ha construido')
    !    
    !n es el numero de celdas activas del sistema. Tambien puede ser interpretado como el numero de variables del sistema de 
    !acuaciones o el numero de filas o columna de la matriz densa de coeficientes.
    !
    n=topo%act
!    n=assert_eq((/SIZE(Teq,1),SIZE(Teq,2),topo%act/),'ensmblT')   !Verifica coherencia en el tamaño de arreglos
    !
    !Propiedades generales del sistema    
    rtact  = trim(prjctrt)//'/sanpck.act'   !Archivo de localizaciones de celdas activas
    rthe   = trim(prjctrt)//'/sanpck.dph'   !Archivo de localizaciones de celdas de altura externa
    rttrns = trim(prjctrt)//'/sanpck.cvt'   !Archivo de conductividad de bloque
    OPEN(ua, FILE=rtact, STATUS='OLD', ACCESS='DIRECT', RECL=4*6, FORM='UNFORMATTED')
    OPEN(uhe, FILE=rthe, STATUS='OLD', ACCESS='DIRECT', RECL=4*6, FORM='UNFORMATTED')
    OPEN(u, FILE=rttrns, STATUS='OLD', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
    !
    !Propiedades de las CC de niveles externos    
    rtl     = trim(prjctrt)//'/sanpck.rvl'   !Archivo de longitudes del rio
    rtb     = trim(prjctrt)//'/sanpck.rvb'   !Archivo de anchos del rio
    rtdpth  = trim(prjctrt)//'/sanpck.rdh'   !Archivo de espesor en el lecho de rio
    rtcndtnc= trim(prjctrt)//'/sanpck.rky'   !Archivo de conductividades de lecho de rio
    OPEN(ul, FILE=rtl, STATUS='OLD', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
    OPEN(ub, FILE=rtb, STATUS='OLD', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
    OPEN(w,  FILE=rtdpth, STATUS='OLD', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
    OPEN(v,  FILE=rtcndtnc, STATUS='OLD', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
    !
    dim=topo%cth
    ALLOCATE(Teq(n,n));                         Teq=0.0_dp
    ALLOCATE(qhext(n));                         qhext=0.0_dp
    ALLOCATE(topheconect(n));                   topheconect=0
    ALLOCATE(ijheconect(dim,6));                ijheconect=0
    ALLOCATE(qhcte(dim,1));                     qhcte=0.0_dp
    !
    diagonal: DO a=1,n                            !Primero se recorre la diagonal principal
        !
        i=topo%ijactv(a,1); j=topo%ijactv(a,2); k=topo%ijactv(a,3) !Leo de la  topologia las posiciones
        p=0; Ka=0.0_dp; Kb=0.0_dp; Kc=0.0_dp; Kd=0.0_dp; Ke=0.0_dp; Kf=0.0_dp; Kg=0.0_dp                 
        CALL gtvlr_bn(u,dmngmtry,i,j,k,Ka)
        CALL gt_delta(dmngmtry,i,j,k,dx,dy,dz)  
        !
        !Flujo que entra sobre la cara X
        IF (i/=1) THEN           !Bloque anterior, siempre que no esta en la primera fila de la malla
            CALL gtvlr_bn(ua,dmngmtry,i-1,j,k,stts)  
            IF( (stts==1).OR.(stts==2) ) THEN
                CALL gtvlr_bn(u,dmngmtry,i-1,j,k,Kb)
                CALL gt_distance(dmngmtry,i,j,k,'i-1',dist)
            END IF
            qub=(dy*dz)/(dist)*armncmn(Ka,Kb)
!            CALL changehcterow(p,stts,qub,nheconect,topheconect,ijheconect,qhcte,a,i,j,1)
            Teq(a,a)=Teq(a,a)-qub
        END IF   
        
        !
        !Flujo que sale de la cara X
        IF ( i/=dmngmtry%rws ) THEN !Bloque superior, bloque superior siempre que no est� en la �ltima fila de la malla
            CALL gtvlr_bn(ua,dmngmtry,i+1,j,k,stts)  
            IF( (stts==1).OR.(stts==2) ) THEN
                CALL gtvlr_bn(u,dmngmtry,i+1,j,k,Kc)
                CALL gt_distance(dmngmtry,i,j,k,'i+1',dist)
            END IF
            qub=(dy*dz)/(dist)*armncmn(Ka,Kc)
!            CALL changehcterow(p,stts,qub,nheconect,topheconect,ijheconect,qhcte,a,i,j,2)
            Teq(a,a)=Teq(a,a)-qub 
        END IF   
        !
        !Flujo que entra a la cara Y
        IF ( j/=1 ) THEN          !Bloque izquierdo, siempre que no est� en la primera columna de la malla
            CALL gtvlr_bn(ua,dmngmtry,i,j-1,k,stts)  
            IF( (stts==1).OR.(stts==2) ) THEN
                CALL gtvlr_bn(u,dmngmtry,i,j-1,k,Kd)
                CALL gt_distance(dmngmtry,i,j,k,'j-1',dist)
            END IF
            qub=(dx*dz)/(dist)*armncmn(Ka,Kd)
!            CALL changehcterow(p,stts,qub,nheconect,topheconect,ijheconect,qhcte,a,i,j,3)
            Teq(a,a)=Teq(a,a)-qub 
        END IF   
        !
        !Flujo que sale a la cara Y
        IF (j/=dmngmtry%clmns) THEN !Bloque derecho, siempre que no est� en la �ltima columna de la malla
            CALL gtvlr_bn(ua,dmngmtry,i,j+1,k,stts)  
            IF( (stts==1).OR.(stts==2) ) THEN
                CALL gtvlr_bn(u,dmngmtry,i,j+1,k,Ke)
                CALL gt_distance(dmngmtry,i,j,k,'j+1',dist)
            END IF
            qub=(dx*dz)/(dist)*armncmn(Ka,Ke)
!            CALL changehcterow(p,stts,qub,nheconect,topheconect,ijheconect,qhcte,a,i,j,4)
            Teq(a,a)=Teq(a,a)-qub 
        END IF 
        !
        !Flujo que entra a la cara Z
        IF (k/=1) THEN !Bloque niferior, siempre que no est� en la ultima columna de la malla
            CALL gtvlr_bn(ua,dmngmtry,i,j,k-1,stts)  
            IF( (stts==1).OR.(stts==2) ) THEN
                CALL gtvlr_bn(u,dmngmtry,i,j,k-1,Kf)
                CALL gt_distance(dmngmtry,i,j,k,'k-1',dist)
            END IF
            qub=(dy*dx)/(dist)*armncmn(Ka,Kf)
!            CALL changehcterow(p,stts,qub,nheconect,topheconect,ijheconect,qhcte,a,i,j,5)
            Teq(a,a)=Teq(a,a)-qub 
        END IF 
        !
        !Fujo que sale de la cara Z
        IF (k/=dmngmtry%lvls) THEN !Bloque superior, siempre que no est� en la �ltima columna de la malla
            CALL gtvlr_bn(ua,dmngmtry,i,j,k+1,stts)  
            IF( (stts==1).OR.(stts==2) ) THEN
                CALL gtvlr_bn(u,dmngmtry,i,j,k+1,Kg)
                CALL gt_distance(dmngmtry,i,j,k,'k+1',dist)
            END IF
            qub=(dy*dx)/(dist)*armncmn(Ka,Kg)
!            CALL changehcterow(p,stts,qub,nheconect,topheconect,ijheconect,qhcte,a,i,j,6)
            Teq(a,a)=Teq(a,a)-qub 
        END IF 
        !
        !Flujo asociado a la existencia de una CC tipo 3 o una altura externa.
        !DUDA: ¡Por que interesa almacenar qhext?
        CALL gtvlr_bn(ua,dmngmtry,i,j,k,stts)
        CALL gtvlr_bn(uhe,dmngmtry,i,j,k,sthe)
        conection: IF ( (sthe==3).AND.(stts==1) ) THEN
            CALL gtvlr_bn(w,dmngmtry,i,j,k,E)
            CALL gtvlr_bn(v,dmngmtry,i,j,k,Kr)
            CALL gtvlr_bn(ul,dmngmtry,i,j,k,Le)
            CALL gtvlr_bn(ub,dmngmtry,i,j,k,Wd)
            qub=Kr*Le*Wd/E; qhext(a)=qub !Componente externa de flujo ejercida por la he no nula en el bloque
            Teq(a,a)=Teq(a,a)-qub
        END IF conection
        ! 
    END DO diagonal
    !
    !Ahora se recorre la triangular superior y se utiliza la simetr�a
    activeout: DO a=1,n
        activein: DO b=a+1,n
            
            !Leo de la  topologia las posiciones
            i=topo%ijactv(a,1); j=topo%ijactv(a,2); k=topo%ijactv(a,3)
            l=topo%ijactv(b,1); m=topo%ijactv(b,2); o=topo%ijactv(b,3)
            !Calculo las dimensiones del bloque i,j,k
            CALL gt_delta(dmngmtry,i,j,k,dx,dy,dz)  
            !Conductividad de bloque central
            CALL gtvlr_bn(u,dmngmtry,i,j,k,Ka)
            
            IF ((j==m).AND.(k==o).AND.((i==1+l).OR.(i==l-1))) THEN
            !Transmisividad equivalente para bloques adyacentes activos en una misma columna    
                CALL gtvlr_bn(u,dmngmtry,l,m,o,Kb)
                !
                IF (i==l+1) THEN
                    CALL gt_distance(dmngmtry,i,j,k,'i-1',dist)
                ELSE IF(i==l-1) THEN
                    CALL gt_distance(dmngmtry,i,j,k,'i+1',dist)
                END IF 
                Teq(a,b)=(dy*dz)/dist*armncmn(Ka,Kb)
                !
            ELSE IF ((i==l).AND.(k==o).AND.((j==m+l).OR.(j==m-1))) THEN
            !Transmisividad equivalente para bloques adyacentes activos en una misma columna    
                CALL gtvlr_bn(u,dmngmtry,l,m,o,Kc)
                !
                IF (j==m+1) THEN
                    CALL gt_distance(dmngmtry,i,j,k,'j-1',dist)
                ELSE IF(j==m-1) THEN
                    CALL gt_distance(dmngmtry,i,j,k,'j+1',dist)
                END IF 
                Teq(a,b)=(dx*dz)/dist*armncmn(Ka,Kc)
                !
            ELSE IF ((i==l).AND.(j==m).AND.((k==o+l).OR.(k==o-1))) THEN
            !Transmisividad equivalente para bloques adyacentes activos en una misma columna    
                CALL gtvlr_bn(u,dmngmtry,l,m,o,Kd)
                !
                IF (k==o+1) THEN
                    CALL gt_distance(dmngmtry,i,j,k,'k-1',dist)
                ELSE IF(k==o-1) THEN
                    CALL gt_distance(dmngmtry,i,j,k,'k+1',dist)
                END IF 
                Teq(a,b)=(dx*dy)/dist*armncmn(Ka,Kd)
                
            END IF
            !            
            Teq(b,a)=Teq(a,b)   !Por simetria
            !
        END DO activein        
    END DO activeout
    !
    CLOSE(u); CLOSE(ua); CLOSE(v); CLOSE(uhe)
    CLOSE(w); CLOSE(ul); CLOSE(ub)
    !
    WRITE (*,*) '*---|Finaliza construccion matriz densa de conductividades|---*'
    !
END SUBROUTINE isoTdns_bn
!BL
! Rutina con matrices de Petsc basado en isoTdns_bn, sería ideal que fuese basada en isoTspr_bn cunando está funcione.
SUBROUTINE isoTdns_bn_petsc(prjctrt,topo,dmngmtry,Teq,topheconect,ijheconect,qhcte,qhext)
    USE anisotype
    USE geometry
    USE iodata,     ONLY : gtvlr_bn
    USE isosystem,  ONLY : isochangehcterow
    USE util
    !
    IMPLICIT NONE
    !
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscvec.h>
    !
    PetscErrorCode                                                 :: ierr
    PetscMPIInt                                                    :: tprocess, iprocess,status(MPI_STATUS_SIZE)
    Mat, INTENT(OUT)                                               :: Teq
    !
    CHARACTER(*), INTENT(IN)                                       :: prjctrt
    TYPE(topology), INTENT(IN)                                     :: topo
    TYPE(domain), INTENT(IN)                                       :: dmngmtry
    INTEGER(I4B), DIMENSION(:),ALLOCATABLE, INTENT(OUT)            :: topheconect
    INTEGER(I4B), DIMENSION(:,:),ALLOCATABLE, INTENT(OUT)          :: ijheconect
    REAL(DP), DIMENSION(:,:),ALLOCATABLE, INTENT(OUT)              :: qhcte
    REAL(DP), DIMENSION(:),ALLOCATABLE, INTENT(OUT)                :: qhext
    !
    REAL(DP)       :: Ka,Kb,Kc,Kd,Ke,Kf,Kg,Kr,E,Le,Wd,dx,dy,dz,qub,dist
    CHARACTER(200) :: rtact,rthe,rttrns,rtcndtnc,rtdpth,rtl,rtb
    INTEGER(I4B)   :: a,b,i,j,k,m,l,o,u,v,w,ua,ul,ub,uhe,stts,sthe,n,nheconect,p,dim
    !
    PARAMETER(ua=37);  PARAMETER(u=65)
    PARAMETER(v=58);   PARAMETER(w=33)
    PARAMETER(ul=38);  PARAMETER(ub=39)
    PARAMETER(uhe=74)
    nheconect=0; p=0; 
    !
    CALL MPI_comm_rank(MPI_COMM_WORLD,iprocess,ierr)
    CALL MPI_comm_size(MPI_COMM_WORLD,tprocess,ierr)
    CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
    !
    IF (iprocess==0) THEN
        WRITE (*,*) '*---|Inicia construccion matriz densa de conductividades|---*'
    END IF
    !
    !Revisa si la topologia y geometria se han construido en todos los procesadores
    IF(dmngmtry%exist.EQV..false.)CALL error('En ensmblT: La geometria no se ha construido')
    IF(topo%hastopo.EQV..false.)CALL error('En ensmblT: La topologia no se ha construido')
    !    
    !n es el numero de celdas activas del sistema. Tambien puede ser interpretado como el numero de variables del sistema de 
    !acuaciones o el numero de filas o columna de la matriz densa de coeficientes.
    !
    n=topo%act
    ! Se crea la matriz de Petsc
    CALL MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,n,PETSC_DETERMINE,PETSC_NULL_INTEGER,PETSC_DETERMINE,PETSC_NULL_INTEGER,Teq,ierr)
    CALL MatSetOption(Teq,MAT_SYMMETRIC,PETSC_TRUE,ierr)
    ! Este comando permite dar opciones a la matriz desde la linea de comandos, para ver estas opciones ejecutar el programa adicinando "-help" al final
    ! CALL MatSetFromOptions(Teq,ierr)

!    n=assert_eq((/SIZE(Teq,1),SIZE(Teq,2),topo%act/),'ensmblT')   !Verifica coherencia en el tamaño de arreglos
    !
    IF (iprocess==0) THEN
        !Propiedades generales del sistema    
        rtact  = trim(prjctrt)//'/sanpck.act'   !Archivo de localizaciones de celdas activas
        rthe   = trim(prjctrt)//'/sanpck.dph'   !Archivo de localizaciones de celdas de altura externa
        rttrns = trim(prjctrt)//'/sanpck.cvt'   !Archivo de conductividad de bloque
        OPEN(ua, FILE=rtact, STATUS='OLD', ACCESS='DIRECT', RECL=4*6, FORM='UNFORMATTED')
        OPEN(uhe, FILE=rthe, STATUS='OLD', ACCESS='DIRECT', RECL=4*6, FORM='UNFORMATTED')
        OPEN(u, FILE=rttrns, STATUS='OLD', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
        !
        !Propiedades de las CC de niveles externos    
        rtl     = trim(prjctrt)//'/sanpck.rvl'   !Archivo de longitudes del rio
        rtb     = trim(prjctrt)//'/sanpck.rvb'   !Archivo de anchos del rio
        rtdpth  = trim(prjctrt)//'/sanpck.rdh'   !Archivo de espesor en el lecho de rio
        rtcndtnc= trim(prjctrt)//'/sanpck.rky'   !Archivo de conductividades de lecho de rio
        OPEN(ul, FILE=rtl, STATUS='OLD', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
        OPEN(ub, FILE=rtb, STATUS='OLD', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
        OPEN(w,  FILE=rtdpth, STATUS='OLD', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
        OPEN(v,  FILE=rtcndtnc, STATUS='OLD', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
        !
        dim=topo%cth
        ! ALLOCATE(Teq(n,n));                         Teq=0.0_dp
        ALLOCATE(qhext(n));                         qhext=0.0_dp
        ALLOCATE(topheconect(n));                   topheconect=0
        ALLOCATE(ijheconect(dim,6));                ijheconect=0
        ALLOCATE(qhcte(dim,1));                     qhcte=0.0_dp
        !
        diagonal: DO a=1,n                            !Primero se recorre la diagonal principal
            !
            i=topo%ijactv(a,1); j=topo%ijactv(a,2); k=topo%ijactv(a,3) !Leo de la  topologia las posiciones
            p=0; Ka=0.0_dp; Kb=0.0_dp; Kc=0.0_dp; Kd=0.0_dp; Ke=0.0_dp; Kf=0.0_dp; Kg=0.0_dp                 
            CALL gtvlr_bn(u,dmngmtry,i,j,k,Ka)
            CALL gt_delta(dmngmtry,i,j,k,dx,dy,dz)  
            !
            !Flujo que entra sobre la cara X
            IF (i/=1) THEN           !Bloque anterior, siempre que no esta en la primera fila de la malla
                CALL gtvlr_bn(ua,dmngmtry,i-1,j,k,stts)  
                IF( (stts==1).OR.(stts==2) ) THEN
                    CALL gtvlr_bn(u,dmngmtry,i-1,j,k,Kb)
                    CALL gt_distance(dmngmtry,i,j,k,'i-1',dist)
                END IF
                qub=(dy*dz)/(dist)*armncmn(Ka,Kb)
    !            CALL changehcterow(p,stts,qub,nheconect,topheconect,ijheconect,qhcte,a,i,j,1)
                ! Teq(a,a)=Teq(a,a)-qub
                CALL MatSetValues(Teq,1,a-1,1,a-1,-qub,ADD_VALUES,ierr)
            END IF   
            
            !
            !Flujo que sale de la cara X
            IF ( i/=dmngmtry%rws ) THEN !Bloque superior, bloque superior siempre que no est� en la �ltima fila de la malla
                CALL gtvlr_bn(ua,dmngmtry,i+1,j,k,stts)  
                IF( (stts==1).OR.(stts==2) ) THEN
                    CALL gtvlr_bn(u,dmngmtry,i+1,j,k,Kc)
                    CALL gt_distance(dmngmtry,i,j,k,'i+1',dist)
                END IF
                qub=(dy*dz)/(dist)*armncmn(Ka,Kc)
    !            CALL changehcterow(p,stts,qub,nheconect,topheconect,ijheconect,qhcte,a,i,j,2)
                ! Teq(a,a)=Teq(a,a)-qub 
                CALL MatSetValues(Teq,1,a-1,1,a-1,-qub,ADD_VALUES,ierr)
            END IF   
            !
            !Flujo que entra a la cara Y
            IF ( j/=1 ) THEN          !Bloque izquierdo, siempre que no est� en la primera columna de la malla
                CALL gtvlr_bn(ua,dmngmtry,i,j-1,k,stts)  
                IF( (stts==1).OR.(stts==2) ) THEN
                    CALL gtvlr_bn(u,dmngmtry,i,j-1,k,Kd)
                    CALL gt_distance(dmngmtry,i,j,k,'j-1',dist)
                END IF
                qub=(dx*dz)/(dist)*armncmn(Ka,Kd)
    !            CALL changehcterow(p,stts,qub,nheconect,topheconect,ijheconect,qhcte,a,i,j,3)
                ! Teq(a,a)=Teq(a,a)-qub 
                CALL MatSetValues(Teq,1,a-1,1,a-1,-qub,ADD_VALUES,ierr)
            END IF   
            !
            !Flujo que sale a la cara Y
            IF (j/=dmngmtry%clmns) THEN !Bloque derecho, siempre que no est� en la �ltima columna de la malla
                CALL gtvlr_bn(ua,dmngmtry,i,j+1,k,stts)  
                IF( (stts==1).OR.(stts==2) ) THEN
                    CALL gtvlr_bn(u,dmngmtry,i,j+1,k,Ke)
                    CALL gt_distance(dmngmtry,i,j,k,'j+1',dist)
                END IF
                qub=(dx*dz)/(dist)*armncmn(Ka,Ke)
    !            CALL changehcterow(p,stts,qub,nheconect,topheconect,ijheconect,qhcte,a,i,j,4)
                ! Teq(a,a)=Teq(a,a)-qub 
                CALL MatSetValues(Teq,1,a-1,1,a-1,-qub,ADD_VALUES,ierr)
            END IF 
            !
            !Flujo que entra a la cara Z
            IF (k/=1) THEN !Bloque niferior, siempre que no est� en la ultima columna de la malla
                CALL gtvlr_bn(ua,dmngmtry,i,j,k-1,stts)  
                IF( (stts==1).OR.(stts==2) ) THEN
                    CALL gtvlr_bn(u,dmngmtry,i,j,k-1,Kf)
                    CALL gt_distance(dmngmtry,i,j,k,'k-1',dist)
                END IF
                qub=(dy*dx)/(dist)*armncmn(Ka,Kf)
    !            CALL changehcterow(p,stts,qub,nheconect,topheconect,ijheconect,qhcte,a,i,j,5)
                ! Teq(a,a)=Teq(a,a)-qub 
                CALL MatSetValues(Teq,1,a-1,1,a-1,-qub,ADD_VALUES,ierr)
            END IF 
            !
            !Fujo que sale de la cara Z
            IF (k/=dmngmtry%lvls) THEN !Bloque superior, siempre que no est� en la �ltima columna de la malla
                CALL gtvlr_bn(ua,dmngmtry,i,j,k+1,stts)  
                IF( (stts==1).OR.(stts==2) ) THEN
                    CALL gtvlr_bn(u,dmngmtry,i,j,k+1,Kg)
                    CALL gt_distance(dmngmtry,i,j,k,'k+1',dist)
                END IF
                qub=(dy*dx)/(dist)*armncmn(Ka,Kg)
    !            CALL changehcterow(p,stts,qub,nheconect,topheconect,ijheconect,qhcte,a,i,j,6)
                ! Teq(a,a)=Teq(a,a)-qub 
                CALL MatSetValues(Teq,1,a-1,1,a-1,-qub,ADD_VALUES,ierr)
            END IF 
            !
            !Flujo asociado a la existencia de una CC tipo 3 o una altura externa.
            !DUDA: ¡Por que interesa almacenar qhext?
            CALL gtvlr_bn(ua,dmngmtry,i,j,k,stts)
            CALL gtvlr_bn(uhe,dmngmtry,i,j,k,sthe)
            conection: IF ( (sthe==3).AND.(stts==1) ) THEN
                CALL gtvlr_bn(w,dmngmtry,i,j,k,E)
                CALL gtvlr_bn(v,dmngmtry,i,j,k,Kr)
                CALL gtvlr_bn(ul,dmngmtry,i,j,k,Le)
                CALL gtvlr_bn(ub,dmngmtry,i,j,k,Wd)
                qub=Kr*Le*Wd/E; qhext(a)=qub !Componente externa de flujo ejercida por la he no nula en el bloque
                ! Teq(a,a)=Teq(a,a)-qub
                CALL MatSetValues(Teq,1,a-1,1,a-1,-qub,ADD_VALUES,ierr)
            END IF conection
            ! 
        END DO diagonal
        !
        !Ahora se recorre la triangular superior y se utiliza la simetr�a
        activeout: DO a=1,n
            activein: DO b=a+1,n
                
                !Leo de la  topologia las posiciones
                i=topo%ijactv(a,1); j=topo%ijactv(a,2); k=topo%ijactv(a,3)
                l=topo%ijactv(b,1); m=topo%ijactv(b,2); o=topo%ijactv(b,3)
                !Calculo las dimensiones del bloque i,j,k
                CALL gt_delta(dmngmtry,i,j,k,dx,dy,dz)  
                !Conductividad de bloque central
                CALL gtvlr_bn(u,dmngmtry,i,j,k,Ka)
                
                IF ((j==m).AND.(k==o).AND.((i==1+l).OR.(i==l-1))) THEN
                !Transmisividad equivalente para bloques adyacentes activos en una misma columna    
                    CALL gtvlr_bn(u,dmngmtry,l,m,o,Kb)
                    !
                    IF (i==l+1) THEN
                        CALL gt_distance(dmngmtry,i,j,k,'i-1',dist)
                    ELSE IF(i==l-1) THEN
                        CALL gt_distance(dmngmtry,i,j,k,'i+1',dist)
                    END IF 
                    ! Teq(a,b)=(dy*dz)/dist*armncmn(Ka,Kb)
                    CALL MatSetValues(Teq,1,a-1,1,b-1,(dy*dz)/dist*armncmn(Ka,Kb),ADD_VALUES,ierr)
                    !
                ELSE IF ((i==l).AND.(k==o).AND.((j==m+l).OR.(j==m-1))) THEN
                !Transmisividad equivalente para bloques adyacentes activos en una misma columna    
                    CALL gtvlr_bn(u,dmngmtry,l,m,o,Kc)
                    !
                    IF (j==m+1) THEN
                        CALL gt_distance(dmngmtry,i,j,k,'j-1',dist)
                    ELSE IF(j==m-1) THEN
                        CALL gt_distance(dmngmtry,i,j,k,'j+1',dist)
                    END IF 
                    ! Teq(a,b)=(dx*dz)/dist*armncmn(Ka,Kc)
                    CALL MatSetValues(Teq,1,a-1,1,b-1,(dx*dz)/dist*armncmn(Ka,Kc),ADD_VALUES,ierr)
                    !
                ELSE IF ((i==l).AND.(j==m).AND.((k==o+l).OR.(k==o-1))) THEN
                !Transmisividad equivalente para bloques adyacentes activos en una misma columna    
                    CALL gtvlr_bn(u,dmngmtry,l,m,o,Kd)
                    !
                    IF (k==o+1) THEN
                        CALL gt_distance(dmngmtry,i,j,k,'k-1',dist)
                    ELSE IF(k==o-1) THEN
                        CALL gt_distance(dmngmtry,i,j,k,'k+1',dist)
                    END IF 
                    ! Teq(a,b)=(dx*dy)/dist*armncmn(Ka,Kd)
                    CALL MatSetValues(Teq,1,a-1,1,b-1,(dx*dy)/dist*armncmn(Ka,Kd),ADD_VALUES,ierr)
                    
                END IF
                !            
                ! Teq(b,a)=Teq(a,b)   !Por simetria
                !
            END DO activein        
        END DO activeout
        !
        CLOSE(u); CLOSE(ua); CLOSE(v); CLOSE(uhe)
        CLOSE(w); CLOSE(ul); CLOSE(ub)
        !
    END IF
    !
    CALL MatAssemblyBegin(Teq,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyEnd(Teq,MAT_FINAL_ASSEMBLY,ierr)
    !
    !
    CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
    IF (iprocess==0) THEN
        WRITE (*,*) '*---|Finaliza construccion matriz densa de conductividades|---*'
    END IF
    !
END SUBROUTINE isoTdns_bn_petsc
!BL
SUBROUTINE isoTsprc_bn(prjctrt,topo,dmngmtry,steq,topheconect,ijheconect,qhcte,qhext)
    USE anisotype
    USE geometry
    USE iodata,     ONLY : gtvlr_bn
    USE isosystem,  ONLY : isochangehcterow
    USE util
    !
    IMPLICIT NONE
    !
    CHARACTER(*), INTENT(IN)                                        :: prjctrt
    TYPE(topology), INTENT(IN)                                      :: topo
    TYPE(domain), INTENT(IN)                                        :: dmngmtry
    TYPE(dispersed), INTENT(OUT)                                    :: steq
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE, INTENT(OUT)            :: topheconect
    INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT)          :: ijheconect
    REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT)              :: qhcte
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT)                :: qhext
    !
    REAL(DP), DIMENSION(:), ALLOCATABLE, TARGET     :: val
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE, TARGET :: irow,jcol
    REAL(DP)                                        :: Ka,Kb,Kc,Kd,Ke,Kf,Kg,Kr,Esp,Le,Wd,dx,dy,dz,dist
    REAL(DP)                                        :: qab,qac,qad,qae,qaf,qag,qhe
    CHARACTER(200)                                  :: rtact,rthc,rthe,rttrns,rtcndtnc,rtdpth,rtl,rtb,rttop,rttoc,rttoe
    INTEGER(I4B)                                    :: u,ua,uhc,uhe,ul,ub,v,w,utop,utoe,utoc
    INTEGER(I4B)                                    :: a,b,c,d,e,f,g,i,j,k,m,l,z,stts,sthe,n,nheconect,p,cont,dim
    !
    PARAMETER(ua=37,u=65,v=58,w=33,ul=38,ub=39,uhe=74,utop=21,utoe=22,utoc=23)
    nheconect=0; p=0;
    !
    WRITE (*,*) '*---|Inicia construccion matriz dispersa de conductividades|---*'
    !
    !Revisa si la topologia se ha construido previamente. De no estarlo
    !se arroja un error y el programa se detiene
    IF(dmngmtry%exist.EQV..false.)CALL error('En ensmblT: La geometria no se ha construido')
    IF(topo%hastopo.EQV..false.)CALL error('En ensmblT: La topologia no se ha construido')
    !
    n=topo%act; cont=0              
    steq%n=n                        !Atributo de variable tipo de indica numero de celdas activas
    steq%len=7*n                    !Tamaño maximo de la matriz dispersa. Por cada celda activa hay maximo seis componentes de flujo
    ALLOCATE(steq%val(7*n),steq%irow(7*n),steq%jcol(7*n))
    steq%val=0.0_dp; steq%irow=0; steq%jcol=0
    !
    !Propiedades generales del acuifero    
    rtact  = trim(prjctrt)//'/sanpck.act'   !Archivo de localizaciones de celdas activas
    rthe   = trim(prjctrt)//'/sanpck.dph'   !Archivo de localizaciones de celdas de altura piezometrica exterior    
    rttop  = trim(prjctrt)//'/sanpck.top'   !Archivo que contiene indices topologicos para bloques activos
    rttoc  = trim(prjctrt)//'/sanpck.toc'   !Archivo que contiene indices topologicos para CC de nivel impuesto
    rttoe  = trim(prjctrt)//'/sanpck.toe'   !Archivo que contiene indices topologicos para CC de nivel externo
    rttrns = trim(prjctrt)//'/sanpck.cvt'   !Archivo de transmisividad de los bloques
    OPEN(ua, FILE=rtact, STATUS='OLD', ACCESS='DIRECT', RECL=4*6, FORM='UNFORMATTED')
    OPEN(uhe, FILE=rthe, STATUS='OLD', ACCESS='DIRECT', RECL=4*6, FORM='UNFORMATTED')
    OPEN(utop, FILE=rttop, STATUS='OLD', ACCESS='DIRECT', RECL=4*10, FORM='UNFORMATTED')
    OPEN(utoc, FILE=rttoc, STATUS='OLD', ACCESS='DIRECT', RECL=4*10, FORM='UNFORMATTED')
    OPEN(utoe, FILE=rttoe, STATUS='OLD', ACCESS='DIRECT', RECL=4*10, FORM='UNFORMATTED')
    OPEN(u, FILE=rttrns, STATUS='OLD', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
    !
    !Propiedades de las CC de niveles externos    
    rtl     = trim(prjctrt)//'/sanpck.rvl'   !Archivo de longitudes del rio
    rtb     = trim(prjctrt)//'/sanpck.rvb'   !Archivo de anchos del rio
    rtdpth  = trim(prjctrt)//'/sanpck.rdh'   !Archivo de espesor en el lecho de rio
    rtcndtnc= trim(prjctrt)//'/sanpck.rky'   !Archivo de conductividades de lecho de rio
    OPEN(ul, FILE=rtl, STATUS='OLD', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
    OPEN(ub, FILE=rtb, STATUS='OLD', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
    OPEN(w,  FILE=rtdpth, STATUS='OLD', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
    OPEN(v,  FILE=rtcndtnc, STATUS='OLD', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
    !
    dim=SIZE(topo%ijcth,1)
    ALLOCATE(qhext(n));                         qhext=0.0_dp
    ALLOCATE(topheconect(n));                   topheconect=0
    ALLOCATE(ijheconect(dim,6));                ijheconect=0
    ALLOCATE(qhcte(dim,1));                     qhcte=0.0_dp
    !
    diagonal: DO a=1,n                            
        !
        i=topo%ijactv(a,1); j=topo%ijactv(a,2); k=topo%ijactv(a,3)  !Leo de la  topologia las posiciones
        qab=0.0_dp; qac=0.0_dp; qad=0.0_dp; qae=0.0_dp; qaf=0.0_dp; qag=0.0_dp
        Ka=0.0_dp; Kb=0.0_dp; Kc=0.0_dp; Kd=0.0_dp; Ke=0.0_dp; Kf=0.0_dp; Kg=0.0_dp
        p=0; b=0; c=0; d=0; e=0; f=0; g=0                
        CALL gtvlr_bn(u,dmngmtry,i,j,k,Ka) !Calculo conductividad hidraulica del bloque central
        CALL gt_delta(dmngmtry,i,j,k,dx,dy,dz) !Calculo longitud de los lados de bloque central
        !
        IF (i/=1) THEN           !Bloque anterior en X, siempre que no este en la primera columna de la malla
            CALL gtvlr_bn(ua,dmngmtry,i-1,j,k,stts) 
            IF( (stts==1).OR.(stts==2) ) THEN
                CALL gtvlr_bn(u,dmngmtry,i-1,j,k,Kb)
                CALL gt_distance(dmngmtry,i,j,k,'i-1',dist)
            END IF
            IF (stts==1) THEN
                CALL gtvlr_bn(utop,dmngmtry,i-1,j,k,b)
            END IF
            qab=(dy*dz)/(dist)*armncmn(Ka,Kb)
!            CALL changehcterow(p,stts,qab,nheconect,topheconect,ijheconect,qhcte,a,i,j,1)            
        END IF
        !
        IF ( i/=dmngmtry%clmns ) THEN !Bloque superior en X, bloque superior siempre que no este en la ultimacolumna de la malla
            CALL gtvlr_bn(ua,dmngmtry,i+1,j,k,stts) 
            IF( (stts==1).OR.(stts==2) ) THEN
                CALL gtvlr_bn(u,dmngmtry,i+1,j,k,Kc)
                CALL gt_distance(dmngmtry,i,j,k,'i+1',dist)
            END IF
            IF (stts==1) THEN
                CALL gtvlr_bn(utop,dmngmtry,i+1,j,k,c)
            END IF
            qac=(dy*dz)/(dist)*armncmn(Ka,Kc)
!            CALL changehcterow(p,stts,qac,nheconect,topheconect,ijheconect,qhcte,a,i,j,2)            
        END IF
        !
        IF ( j/=1 ) THEN          !Bloque anterior en Y, siempre que no este en la primera fila de la malla
            CALL gtvlr_bn(ua,dmngmtry,i,j-1,k,stts)  
            IF( (stts==1).OR.(stts==2) ) THEN  
                CALL gtvlr_bn(u,dmngmtry,i,j-1,k,Kd)
                CALL gt_distance(dmngmtry,i,j,k,'j-1',dist)
            END IF
            IF (stts==1) THEN
                CALL gtvlr_bn(utop,dmngmtry,i,j-1,k,d)
            END IF
            qad=(dx*dz)/(dist)*armncmn(Ka,Kd)
!            CALL changehcterow(p,stts,qad,nheconect,topheconect,ijheconect,qhcte,a,i,j,3)            
        END IF   
        !
        IF (j/=dmngmtry%rws) THEN !Bloque superior en Y, siempre que no este en la ultima fila de la malla
            CALL gtvlr_bn(ua,dmngmtry,i,j+1,k,stts)  
            IF( (stts==1).OR.(stts==2) ) THEN
                CALL gtvlr_bn(u,dmngmtry,i,j+1,k,Ke)
                CALL gt_distance(dmngmtry,i,j,k,'j+1',dist)
            END IF
            IF (stts==1) THEN
                CALL gtvlr_bn(utop,dmngmtry,i,j+1,k,e)
            END IF
            qae=(dx*dz)/(dist)*armncmn(Ka,Ke)
!            CALL changehcterow(p,stts,qae,nheconect,topheconect,ijheconect,qhcte,a,i,j,4)            
        END IF   
        !
        IF ( k/=1 ) THEN          !Bloque anterior en Z, siempre que no este en la primera capa de la malla
            CALL gtvlr_bn(ua,dmngmtry,i,j,k-1,stts)  
            IF( (stts==1).OR.(stts==2) ) THEN
                CALL gtvlr_bn(u,dmngmtry,i,j,k-1,Kf)
                CALL gt_distance(dmngmtry,i,j,k,'k-1',dist)
            END IF
            IF (stts==1) THEN
                CALL gtvlr_bn(utop,dmngmtry,i,j,k-1,f)
            END IF
            qaf=(dx*dy)/(dist)*armncmn(Ka,Kf)
!            CALL changehcterow(p,stts,qad,nheconect,topheconect,ijheconect,qhcte,a,i,j,5)            
        END IF 
        !
        IF (k/=dmngmtry%clmns) THEN !Bloque superior en Z, siempre que no este en la ultima capa de la malla
            CALL gtvlr_bn(ua,dmngmtry,i,j,k+1,stts)    
            IF( (stts==1).OR.(stts==2) ) THEN
                CALL gtvlr_bn(u,dmngmtry,i,j,k+1,Kg)
                CALL gt_distance(dmngmtry,i,j,k,'k+1',dist)
            END IF
            IF (stts==1) THEN
                CALL gtvlr_bn(utop,dmngmtry,i,j,k+1,g)
            END IF
            qag=(dx*dy)/(dist)*armncmn(Ka,Kg)
!            CALL changehcterow(p,stts,qae,nheconect,topheconect,ijheconect,qhcte,a,i,j,6)            
        END IF   
        !
        !Primer bloque del recorrido, inferior Z
        IF(f/=0 .AND. Kf/=0.0_dp)THEN
            cont=cont+1
            steq%jcol(cont)=a
            steq%irow(cont)=f
            steq%val(cont)=qaf
        END IF
        !
        !Segundo bloque del recorrido, inferior Y
        IF(d/=0 .AND. Kd/=0.0_dp)THEN
            cont=cont+1
            steq%jcol(cont)=a
            steq%irow(cont)=d
            steq%val(cont)=qad
        END IF
        !
        !Tercer bloque del recorrido, inferior X
        IF(b/=0 .AND. Kb/=0.0_dp)THEN
            cont=cont+1
            steq%jcol(cont)=a
            steq%irow(cont)=b
            steq%val(cont)=qab
        END IF
        !
        !Cuarto bloque central
        cont=cont+1
        steq%jcol(cont)=a
        steq%irow(cont)=a
        steq%val(cont)=steq%val(cont)-(qab+qac+qad+qae+qaf+qag)
        !
        !Se verifica si la celda activa pertenece a la CC de altura exterior
        CALL gtvlr_bn(ua,dmngmtry,i,j,k,stts)
        CALL gtvlr_bn(uhe,dmngmtry,i,j,k,sthe)
        conection: IF ( (sthe==3).AND.(stts==1) ) THEN
            CALL gtvlr_bn(w,dmngmtry,i,j,k,Esp)
            CALL gtvlr_bn(v,dmngmtry,i,j,k,Kr)
            CALL gtvlr_bn(ul,dmngmtry,i,j,k,Le)
            CALL gtvlr_bn(ub,dmngmtry,i,j,k,Wd)
            qhe=Kr*Le*Wd/Esp; qhext(a)=qhe !Componente externa de flujo ejercida por la altura externa en el bloque
            steq%val(cont)=steq%val(cont)-qhe            
        END IF conection 
        !Quinto bloque del recorrido, superior en X
        IF(c/=0 .AND. Kc/=0.0_dp)THEN
            cont=cont+1
            steq%jcol(cont)=a
            steq%irow(cont)=c
            steq%val(cont)=qac
        END IF    
        !
        !Sexto bloque del recorrido, superior en Y
        IF(e/=0 .AND. Ke/=0.0_dp)THEN
            cont=cont+1
            steq%jcol(cont)=a
            steq%irow(cont)=e
            steq%val(cont)=qae
        END IF
        !
        !Septimo bloque del recorrido, superior en Z
        IF(g/=0 .AND. Kg/=0.0_dp)THEN
            cont=cont+1
            steq%jcol(cont)=a
            steq%irow(cont)=g
            steq%val(cont)=qag
        END IF
    END DO diagonal
    !
    !Probar si saliendo del bucle aLLOcate con punteros
    steq%len=cont
    !
    CLOSE(u); CLOSE(ua); CLOSE(v); CLOSE(uhe)
    CLOSE(w); CLOSE(ul); CLOSE(ub)
    CLOSE(utop); CLOSE(utoe); CLOSE(utoc)
    
    WRITE (*,*) '*---|Finaliza construccion matriz dispersa de conductividades|---*'
    !
END SUBROUTINE isoTsprc_bn
!
!***********************************************************************************************************************
!*******SUBRUTINA QUE ENSAMBLA LA MATRIZ DE ALMACENAMIENTOS MODIFICADA**************************************************
!*****************@author �lvarez-Villa O.D.****************************************************************************
!***********************************************************************************************************************
!
!Subrutina para ensamblar la matriz de almacenamientos. Los
!par�metros que recibe son:
!prjctrt        :Ruta del proyecto de trabajo.
!prjnm          :Nombre del proyecto de trabajo.
!topo          :La topologia del acuifero lineal
!dmngmtry          :La geometria del acuifero
!Salidas:
!SF             :El vector de almacenamientos
!
SUBROUTINE isoSF(prjctrt,topo,dmngmtry,SF)
    USE anisotype
    USE geometry
    USE iodata,        ONLY : gtvlr_bn
    USE util
    !
    IMPLICIT NONE
    !    
    CHARACTER(*), INTENT(IN)                :: prjctrt
    TYPE(topology), INTENT(IN)              :: topo
    TYPE(domain), INTENT(IN)                :: dmngmtry
    REAL(DP), DIMENSION(:), INTENT(INOUT)   :: SF
    !
    REAL(DP)       :: sff,dx,dy,dz
    INTEGER(I4B)   :: a,b,i,j,k,u,n
    CHARACTER(200) :: rtstrg    
    !
    PARAMETER(u=65)
    !
    WRITE (*,*) '*---|Inicia construccion vector de almacenamientos|---*'
    !
    !Revisa si la topologia y geometria se han construido
    IF(dmngmtry%exist.EQV..false.)CALL error('En ensmblT: La geometria no se ha construido')
    IF(topo%hastopo.EQV..false.)CALL error('En ensmblT: La topologia no se ha construido')
    !
    n=assert_eq((/SIZE(SF),topo%act/),'anisoSF')   !Verificar tamaño de los arreglos
    !
    rtstrg=prjctrt//'/sanpck.stg'
    OPEN(u, FILE=rtstrg, STATUS='OLD', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')    
    !        
    activeout: DO a=1,n
        i=topo%ijactv(a,1); j=topo%ijactv(a,2); k=topo%ijactv(a,3)  !Leo de la  topologia las posiciones
        CALL gtvlr_bn(u,dmngmtry,i,j,k,sff)
        CALL gt_delta(dmngmtry,i,j,k,dx,dy,dz)                      !Calculo longitud de los lados de bloque central
        
        SF(a)=dx*dy*dz*sff;
    END DO activeout
    !
    WRITE (*,*) '*---|Finaliza construccion vector de almacenamientos|---*'
    !
END SUBROUTINE isoSF

!***********************************************************************************************************************
!*******************SUBRUTINAS AUXILIARES*******************************************************************************
!*****************@author Alvarez-Villa O.D.****************************************************************************
!***********************************************************************************************************************
!
!Subrutina que inserta y cambia las fila de las matrices indicadoras de conexion
!con CC de altura piezom�trica impuesta. los par�metros son:
!p          :Entero que indica si se inserta fila nueva o se modifica la actual.
!stts       :Entero que indica si se trata de una CC de altura impuesta.
!qub        :Transmisividad equivalente de intercambio.                      
!nheconect  :Entero que indica la fila actual o bloque de modificaci�n.
!topheconect:Arreglo de indicadores de celdas activas donde existe conexion.
!ijhcinect  :Arreglo que indican cuales celdas con CC de himpuesta se conectan con las activas
!qhcte      :Arreglo con las transmisividades equivalentes de las celdas activas con conexion.
!a          :Entero que indica la posicion dentro de la topologia de la celda activa con conexion.
!i          :Localizaci�n matricial de la celda activa, fila.
!j          :Localizaci�n matricial de la celda activa, columna.
!pstn       :Indicador de las celdas de conexion: 1:abajo, 2:arriba, 3:izquierda, 4:derecha
!
SUBROUTINE isochangehcterow(p,stts,qub,nheconect,topheconect,ijheconect,qhcte,a,i,j,pstn)
    USE anisotype
    USE util
    !
    IMPLICIT NONE
    !    
    INTEGER(I4B), INTENT(INOUT)                                       :: p
    INTEGER(I4B), INTENT(IN)                                          :: stts
    REAL(DP), INTENT(IN)                                              :: qub
    INTEGER(I4B), INTENT(INOUT)                                       :: nheconect
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)            :: topheconect
    INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)          :: ijheconect
    REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)              :: qhcte
    INTEGER(I4B), INTENT(IN)                                          :: a,i,j,pstn
    !
    IF(stts==2)THEN
        p=p+1
        IF(p==1)THEN
            nheconect=nheconect+1
            topheconect(a)=nheconect           
        END IF
        ijheconect(nheconect,pstn)=1
        qhcte(nheconect,pstn)=qub
    END IF
    !
END SUBROUTINE isochangehcterow