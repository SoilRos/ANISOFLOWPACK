!-----------------------------------------------------------------------------------------------------------------------------------
!Copyright (c) 2013, Alvarez-Villa O.D., Perez Kevin  GOTTA Ingenieria SAS, All Right Reserved.
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
MODULE anisosystem
    INTERFACE anisoensmblT_bn
        SUBROUTINE anisoTdns_bn(prjctrt,topo,dmngmtry,Teq,topheconect,ijheconect,qhcte,qhext,kxx,kyy,kzz,kxy,kxz,kyz)
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
            REAL(DP), DIMENSION(:),ALLOCATABLE, INTENT(OUT)                :: kxx,kyy,kzz,kxy,kxz,kyz
        END SUBROUTINE anisoTdns_bn
    END INTERFACE anisoensmblT_bn    
    !
END MODULE anisosystem

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

SUBROUTINE anisoTdns_bn(prjctrt,topo,dmngmtry,Teq,topheconect,ijheconect,qhcte,qhext,kxx,kyy,kzz,kxy,kxz,kyz)
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
    REAL(DP), DIMENSION(:),ALLOCATABLE, INTENT(OUT)                :: kxx,kyy,kzz,kxy,kxz,kyz
    !
    REAL(DP)       :: Ka,Kb,Kc,Kd,Ke,Kf,Kg,Kr,E,Le,Wd,dx,dy,dz,qub,dist,delta,delte,delti,delto
    CHARACTER(200) :: rtact,rthe,rttrns,rtcndtnc,rtdpth,rtl,rtb,rtkx,rtky,rtkz,rtxy,rtxz,rtyz
    INTEGER(I4B)   :: a,b,i,j,k,m,l,o,stts,sthe,n,nheconect,p,dim,acto,acti
    INTEGER(I4B)   :: ua,ul,ub,uhe,v,w,utop,utoc,utoe,uc,ud,ue,uf,ug,uh,u
    !
    PARAMETER(ua=37);  PARAMETER(ub=65)
    PARAMETER(v=58);   PARAMETER(w=33)
    PARAMETER(ul=38);  PARAMETER(utoc=39)
    PARAMETER(uhe=74); PARAMETER(utop=66)
    PARAMETER(utoe=35); PARAMETER(uc=44)
    PARAMETER(ud=46); PARAMETER(ue=31)
    PARAMETER(uf=69); PARAMETER(ug=78)
    PARAMETER(uh=42); PARAMETER(U=11)
    nheconect=0; p=0; 
    !
    WRITE (*,*) '*---|Inicia construccion matriz densa con tensor anisotropo de conductividades |---*'
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
    rthe   = trim(prjctrt)//'/sanpck.dph'   !Archivo de localizaciones de celdas de altura extern
    rtkx = trim(prjctrt)//'/sanpck.kxx'     !Archivo con conductividades en direccion xx
    rtky = trim(prjctrt)//'/sanpck.kyy'     !Archivo con conductividades en direccion yy
    rtkz = trim(prjctrt)//'/sanpck.kzz'     !Archivo con conductividades en direccion zz
    rtxy= trim(prjctrt)//'/sanpck.kxy'      !Archivo con conductividades en direccion xy
    rtxz = trim(prjctrt)//'/sanpck.kxz'     !Archivo con conductividades en direccion xz
    rtyz = trim(prjctrt)//'/sanpck.kyz'     !Archivo con conductividades en direccion yz
    rttrns = trim(prjctrt)//'/sanpck.cvt'   !Archivo de conductividad de bloque
    OPEN(ua, FILE=rtact, STATUS='OLD', ACCESS='DIRECT', RECL=4*6, FORM='UNFORMATTED')
    OPEN(uhe, FILE=rthe, STATUS='OLD', ACCESS='DIRECT', RECL=4*6, FORM='UNFORMATTED')
    OPEN(uc, FILE=rtkx, STATUS='OLD', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
    OPEN(ud, FILE=rtky, STATUS='OLD', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
    OPEN(ue, FILE=rtkz, STATUS='OLD', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
    OPEN(uf, FILE=rtxy, STATUS='OLD', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
    OPEN(ug, FILE=rtxz, STATUS='OLD', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
    OPEN(uh, FILE=rtyz, STATUS='OLD', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
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
    CLOSE(ua); CLOSE(uhe); CLOSE(utop); CLOSE(utoc); CLOSE(utoe)
    CLOSE(uc); CLOSE(ud); CLOSE(ue); CLOSE(uf); CLOSE(ug)
    CLOSE(uh); CLOSE(ul); CLOSE(ub); CLOSE(w); CLOSE(v)
    !
    dim=topo%cth
    ALLOCATE(Teq(n,n));                         Teq=0.0_dp
    ALLOCATE(qhext(n));                         qhext=0.0_dp
    ALLOCATE(topheconect(n));                   topheconect=0
    ALLOCATE(ijheconect(dim,6));                ijheconect=0
    ALLOCATE(qhcte(dim,1));                     qhcte=0.0_dp
    !
    !Calculo de terminos de la diagonal principal
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
    !Se calculas las componentes de anisotropia de flujo. Se aprovecha la simetria de la matriz
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
            !Conductividad equivalente para bloques adyacentes activos en una misma columna    
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
            !Conductividad equivalente para bloques adyacentes activos en una misma columna    
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
            !Conductividad equivalente para bloques adyacentes activos en una misma columna    
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
    
    !Se calculan las componentes anisotropas de flujo del tensor
    anisorow:  DO a=1,n
        anisocol: DO b=1,n
            !
            delta=0.0_dp; delte=0.0_dp; delti=0.0_dp; delto=0.0_dp
            acto=0; acti=0
            Ka=0.0_dp; Kb=0.0_dp; Kc=0.0_dp; Kd=0.0_dp; Ke=0.0_dp;
            !Leo de la  topologia las posiciones
            i=topo%ijactv(a,1); j=topo%ijactv(a,2); k=topo%ijactv(a,3)
            l=topo%ijactv(b,1); m=topo%ijactv(b,2); o=topo%ijactv(b,3)
            !Calculo las dimensiones del bloque i,j,k
            CALL gt_delta(dmngmtry,i,j,k,dx,dy,dz)  
            !Conductividad de bloque central
            CALL gtvlr_bn(uf,dmngmtry,i,j,k,ka) !Conductividad Kxy
            CALL gtvlr_bn(ug,dmngmtry,i,j,k,kb) !Conductividad Kxz
            CALL gtvlr_bn(uh,dmngmtry,i,j,k,kc) !Conductividad Kyz
            !
            !Para bloque i-1,j-1,k  #1
            IF ((i==l+1).AND.(j==m+1).AND.(k==o)) THEN
                CALL gtvlr_bn(ua,dmngmtry,i-1,j+1,k,acto) 
                CALL gtvlr_bn(ua,dmngmtry,i+1,j-1,k,acti) 
                CALL gtvlr_bn(uf,dmngmtry,i-1,j,k,kd)!Conductividad Kxy
                CALL gtvlr_bn(uf,dmngmtry,i,j-1,k,ke)!Conductividad Kxy
                CALL gt_distance(dmngmtry,i-1,j,k,'j-1',delta) !delta y(i-1,j-1,k)
                j=j+1 !Artificio para utilizar subrutina get_distance de geometria
                CALL gt_distance(dmngmtry,i-1,j,k,'j-1',delte) !delta y(i-1,j+1,k)
                j=j-1
                CALL gt_distance(dmngmtry,i,j-1,k,'i-1',delti) !delta x(i-1,j-1,k)
                i=i+1 !Artificio para utilizar subrutina get_distance de geometria
                CALL gt_distance(dmngmtry,i,j-1,k,'i-1',delto) !delta x(i+1,j-1,k)
                i=i-1
                IF ((acto/=0).AND.(acti/=0)) THEN 
                    Teq(a,b)=-(1/dx*(armncmn(ka,kd)/(delta+delte))+(1/dy*(armncmn(ka,ke)/(delti+delto))))
                ELSE IF (acto/=0) THEN
                    Teq(a,b)=-(1/dx*(armncmn(ka,kd)/(delta+delte)))
                ELSE IF (acti/=0) THEN
                    Teq(a,b)=-(1/dy*(armncmn(ka,ke)/(delti+delto)))
                ELSE
                    Teq(a,b)=0
                END IF
            !Para bloque i+1,j-1,k   #2
            ELSE IF ((i==l-1).AND.(j==m+1).AND.(k==o)) THEN
                CALL gtvlr_bn(ua,dmngmtry,i+1,j+1,k,acto) 
                CALL gtvlr_bn(ua,dmngmtry,i-1,j-1,k,acti) 
                CALL gtvlr_bn(uf,dmngmtry,i+1,j,k,kd)
                CALL gtvlr_bn(uf,dmngmtry,i,j-1,k,ke)
                CALL gt_distance(dmngmtry,i+1,j,k,'j-1',delta) 
                j=j+1 !Artificio para utilizar subrutina get_distance de geometria
                CALL gt_distance(dmngmtry,i+1,j,k,'j-1',delte)
                j=j-1
                CALL gt_distance(dmngmtry,i,j-1,k,'i-1',delti) 
                i=i+1 !Artificio para utilizar subrutina get_distance de geometria
                CALL gt_distance(dmngmtry,i,j-1,k,'i-1',delto) 
                i=i-1
                IF ((acto/=0).AND.(acti/=0)) THEN 
                    Teq(a,b)=(1/dx*(armncmn(ka,kd)/(delta+delte))+(1/dy*(armncmn(ka,ke)/(delti+delto))))
                ELSE IF (acto/=0) THEN
                    Teq(a,b)=1/dx*(armncmn(ka,kd)/(delta+delte))
                ELSE IF (acti/=0) THEN
                    Teq(a,b)=1/dy*(armncmn(ka,ke)/(delti+delto))
                ELSE
                    Teq(a,b)=0
                END IF
            !Para bloque i+1,j+1,k   #3
            ELSE IF ((i==l-1).AND.(j==m-1).AND.(k==o)) THEN
                CALL gtvlr_bn(ua,dmngmtry,i+1,j-1,k,acto) 
                CALL gtvlr_bn(ua,dmngmtry,i-1,j+1,k,acti) 
                CALL gtvlr_bn(uf,dmngmtry,i+1,j,k,kd)
                CALL gtvlr_bn(uf,dmngmtry,i,j+1,k,ke)
                CALL gt_distance(dmngmtry,i+1,j,k,'j-1',delta) 
                j=j+1 !Artificio para utilizar subrutina get_distance de geometria
                CALL gt_distance(dmngmtry,i+1,j,k,'j-1',delte) 
                j=j-1
                CALL gt_distance(dmngmtry,i,j+1,k,'i-1',delti) !delta x(i+1,j+1,k)
                i=i+1 !Artificio para utilizar subrutina get_distance de geometria
                CALL gt_distance(dmngmtry,i,j+1,k,'i-1',delto) 
                i=i-1
                IF ((acto/=0).AND.(acti/=0)) THEN 
                    Teq(a,b)=-(1/dx*(armncmn(ka,kd)/(delta+delte))+(1/dy*(armncmn(ka,ke)/(delti+delto))))
                ELSE IF (acto/=0) THEN
                    Teq(a,b)=-(1/dx*(armncmn(ka,kd)/(delta+delte)))
                ELSE IF (acti/=0) THEN
                    Teq(a,b)=-(1/dy*(armncmn(ka,ke)/(delti+delto)))
                ELSE
                    Teq(a,b)=0
                END IF
            !Para bloque i-1,j+1,k   #4
            ELSE IF ((i==l+1).AND.(j==m-1).AND.(k==o)) THEN
                CALL gtvlr_bn(ua,dmngmtry,i-1,j-1,k,acto) 
                CALL gtvlr_bn(ua,dmngmtry,i+1,j+1,k,acti) 
                CALL gtvlr_bn(uf,dmngmtry,i-1,j,k,kd)
                CALL gtvlr_bn(uf,dmngmtry,i,j+1,k,ke)
                CALL gt_distance(dmngmtry,i-1,j,k,'j-1',delta) 
                j=j+1 !Artificio para utilizar subrutina get_distance de geometria
                CALL gt_distance(dmngmtry,i-1,j,k,'j-1',delte) 
                j=j-1
                CALL gt_distance(dmngmtry,i,j+1,k,'i-1',delti) 
                i=i+1 !Artificio para utilizar subrutina get_distance de geometria
                CALL gt_distance(dmngmtry,i,j+1,k,'i-1',delto) 
                i=i-1
                IF ((acto/=0).AND.(acti/=0)) THEN 
                    Teq(a,b)=(1/dx*(armncmn(ka,kd)/(delta+delte))+(1/dy*(armncmn(ka,ke)/(delti+delto))))
                ELSE IF (acto/=0) THEN
                    Teq(a,b)=1/dx*(armncmn(ka,kd)/(delta+delte))
                ELSE IF (acti/=0) THEN
                    Teq(a,b)=1/dy*(armncmn(ka,ke)/(delti+delto))
                ELSE
                    Teq(a,b)=0
                END IF
            !Para bloque i-1,j,k-1   #5
            ELSE IF ((i==l+1).AND.(j==m).AND.(k==o+1)) THEN
                CALL gtvlr_bn(ua,dmngmtry,i-1,j,k+1,acto) 
                CALL gtvlr_bn(ua,dmngmtry,i+1,j,k-1,acti) 
                CALL gtvlr_bn(ug,dmngmtry,i-1,j,k,kd)
                CALL gtvlr_bn(ug,dmngmtry,i,j,k-1,ke)
                CALL gt_distance(dmngmtry,i-1,j,k,'k-1',delta) !delta z
                k=k+1 !Artificio para utilizar subrutina get_distance de geometria
                CALL gt_distance(dmngmtry,i-1,j,k,'k-1',delte) !delta z
                k=k-1
                CALL gt_distance(dmngmtry,i,j,k-1,'i-1',delti) !delta x
                i=i+1 !Artificio para utilizar subrutina get_distance de geometria
                CALL gt_distance(dmngmtry,i,j,k-1,'i-1',delto) !delta x
                i=i-1
                IF ((acto/=0).AND.(acti/=0)) THEN 
                    Teq(a,b)=-(1/dx*(armncmn(kb,kd)/(delta+delte))+(1/dz*(armncmn(kb,ke)/(delti+delto))))
                ELSE IF (acto/=0) THEN
                    Teq(a,b)=-(1/dx*(armncmn(kb,kd)/(delta+delte)))
                ELSE IF (acti/=0) THEN
                    Teq(a,b)=-(1/dz*(armncmn(kb,ke)/(delti+delto)))
                ELSE
                    Teq(a,b)=0
                END IF
            !Para bloque i-1,j,k+1   #6
            ELSE IF ((i==l+1).AND.(j==m).AND.(k==o-1)) THEN
                CALL gtvlr_bn(ua,dmngmtry,i-1,j,k-1,acto) 
                CALL gtvlr_bn(ua,dmngmtry,i+1,j,k+1,acti) 
                CALL gtvlr_bn(ug,dmngmtry,i-1,j,k,kd)
                CALL gtvlr_bn(ug,dmngmtry,i,j,k+1,ke)
                CALL gt_distance(dmngmtry,i-1,j,k,'k-1',delta) !delta z
                k=k+1 !Artificio para utilizar subrutina get_distance de geometria
                CALL gt_distance(dmngmtry,i-1,j,k,'k-1',delte) !delta z
                k=k-1
                CALL gt_distance(dmngmtry,i,j,k+1,'i-1',delti) !delta x
                i=i+1 !Artificio para utilizar subrutina get_distance de geometria
                CALL gt_distance(dmngmtry,i,j,k+1,'i-1',delto) !delta x
                i=i-1
                IF ((acto/=0).AND.(acti/=0)) THEN 
                    Teq(a,b)=(1/dx*(armncmn(kb,kd)/(delta+delte))+(1/dz*(armncmn(kb,ke)/(delti+delto))))
                ELSE IF (acto/=0) THEN
                    Teq(a,b)=(1/dx*(armncmn(kb,kd)/(delta+delte)))
                ELSE IF (acti/=0) THEN
                    Teq(a,b)=(1/dz*(armncmn(kb,ke)/(delti+delto)))
                ELSE
                    Teq(a,b)=0
                END IF
            !Para bloque i+1,j,k+1   #7
            ELSE IF ((i==l-1).AND.(j==m).AND.(k==o-1)) THEN
                CALL gtvlr_bn(ua,dmngmtry,i+1,j,k-1,acto) 
                CALL gtvlr_bn(ua,dmngmtry,i-1,j,k+1,acti) 
                CALL gtvlr_bn(ug,dmngmtry,i+1,j,k,kd)
                CALL gtvlr_bn(ug,dmngmtry,i,j,k+1,ke)
                CALL gt_distance(dmngmtry,i+1,j,k,'k-1',delta) !delta z
                k=k+1 !Artificio para utilizar subrutina get_distance de geometria
                CALL gt_distance(dmngmtry,i+1,j,k,'k-1',delte) !delta z
                k=k-1
                CALL gt_distance(dmngmtry,i,j,k+1,'i-1',delti) !delta x
                i=i+1 !Artificio para utilizar subrutina get_distance de geometria
                CALL gt_distance(dmngmtry,i,j,k+1,'i-1',delto) !delta x
                i=i-1
                IF ((acto/=0).AND.(acti/=0)) THEN 
                    Teq(a,b)=-(1/dx*(armncmn(kb,kd)/(delta+delte))+(1/dz*(armncmn(kb,ke)/(delti+delto))))
                ELSE IF (acto/=0) THEN
                    Teq(a,b)=-(1/dx*(armncmn(kb,kd)/(delta+delte)))
                ELSE IF (acti/=0) THEN
                    Teq(a,b)=-(1/dz*(armncmn(kb,ke)/(delti+delto)))
                ELSE
                    Teq(a,b)=0
                END IF
            !Para bloque i+1,j,k-1   #8
            ELSE IF ((i==l-1).AND.(j==m).AND.(k==o+1)) THEN
                CALL gtvlr_bn(ua,dmngmtry,i+1,j,k+1,acto) 
                CALL gtvlr_bn(ua,dmngmtry,i-1,j,k-1,acti) 
                CALL gtvlr_bn(ug,dmngmtry,i+1,j,k,kd)
                CALL gtvlr_bn(ug,dmngmtry,i,j,k-1,ke)
                CALL gt_distance(dmngmtry,i+1,j,k,'k-1',delta) !delta z
                k=k+1 !Artificio para utilizar subrutina get_distance de geometria
                CALL gt_distance(dmngmtry,i+1,j,k,'k-1',delte) !delta z
                k=k-1
                CALL gt_distance(dmngmtry,i,j,k-1,'i-1',delti) !delta x
                i=i+1 !Artificio para utilizar subrutina get_distance de geometria
                CALL gt_distance(dmngmtry,i,j,k-1,'i-1',delto) !delta x
                i=i-1
                IF ((acto/=0).AND.(acti/=0)) THEN 
                    Teq(a,b)=(1/dx*(armncmn(kb,kd)/(delta+delte))+(1/dz*(armncmn(kb,ke)/(delti+delto))))
                ELSE IF (acto/=0) THEN
                    Teq(a,b)=(1/dx*(armncmn(kb,kd)/(delta+delte)))
                ELSE IF (acti/=0) THEN
                    Teq(a,b)=(1/dz*(armncmn(kb,ke)/(delti+delto)))
                ELSE
                    Teq(a,b)=0
                END IF
            !Para bloque i,j-1,k-1   #9
            ELSE IF ((i==l).AND.(j==m+1).AND.(k==o+1)) THEN
                CALL gtvlr_bn(ua,dmngmtry,i,j-1,k+1,acto) 
                CALL gtvlr_bn(ua,dmngmtry,i,j+1,k-1,acti) 
                CALL gtvlr_bn(uh,dmngmtry,i,j-1,k,kd)
                CALL gtvlr_bn(uh,dmngmtry,i,j,k-1,ke)
                CALL gt_distance(dmngmtry,i,j-1,k,'k-1',delta) !delta z
                k=k+1 !Artificio para utilizar subrutina get_distance de geometria
                CALL gt_distance(dmngmtry,i,j-1,k,'k-1',delte) !delta z
                k=k-1
                CALL gt_distance(dmngmtry,i,j,k-1,'j-1',delti) !delta y
                j=j+1 !Artificio para utilizar subrutina get_distance de geometria
                CALL gt_distance(dmngmtry,i,j,k-1,'j-1',delto) !delta y
                j=j-1
                IF ((acto/=0).AND.(acti/=0)) THEN 
                    Teq(a,b)=-(1/dy*(armncmn(kc,kd)/(delta+delte))+(1/dz*(armncmn(kc,ke)/(delti+delto))))
                ELSE IF (acto/=0) THEN
                    Teq(a,b)=-(1/dy*(armncmn(kc,kd)/(delta+delte)))
                ELSE IF (acti/=0) THEN
                    Teq(a,b)=-(1/dz*(armncmn(kc,ke)/(delti+delto)))
                ELSE
                    Teq(a,b)=0
                END IF         
            !Para bloque i,j-1,k+1   #10
            ELSE IF ((i==l).AND.(j==m+1).AND.(k==o-1)) THEN
                CALL gtvlr_bn(ua,dmngmtry,i,j-1,k-1,acto) 
                CALL gtvlr_bn(ua,dmngmtry,i,j+1,k+1,acti) 
                CALL gtvlr_bn(uh,dmngmtry,i,j-1,k,kd)
                CALL gtvlr_bn(uh,dmngmtry,i,j,k+1,ke)
                CALL gt_distance(dmngmtry,i,j-1,k,'k-1',delta) !delta z
                k=k+1 !Artificio para utilizar subrutina get_distance de geometria
                CALL gt_distance(dmngmtry,i,j-1,k,'k-1',delte) !delta z
                k=k-1
                CALL gt_distance(dmngmtry,i,j,k+1,'j-1',delti) !delta y
                j=j+1 !Artificio para utilizar subrutina get_distance de geometria
                CALL gt_distance(dmngmtry,i,j,k+1,'j-1',delto) !delta y
                j=j-1
                IF ((acto/=0).AND.(acti/=0)) THEN 
                    Teq(a,b)=(1/dy*(armncmn(kc,kd)/(delta+delte))+(1/dz*(armncmn(kc,ke)/(delti+delto))))
                ELSE IF (acto/=0) THEN
                    Teq(a,b)=(1/dy*(armncmn(kc,kd)/(delta+delte)))
                ELSE IF (acti/=0) THEN
                    Teq(a,b)=(1/dz*(armncmn(kc,ke)/(delti+delto)))
                ELSE
                    Teq(a,b)=0
                END IF  
            !Para bloque i,j+1,k+1   #11
            ELSE IF ((i==l).AND.(j==m-1).AND.(k==o-1)) THEN
                CALL gtvlr_bn(ua,dmngmtry,i,j+1,k-1,acto) 
                CALL gtvlr_bn(ua,dmngmtry,i,j-1,k+1,acti) 
                CALL gtvlr_bn(uh,dmngmtry,i,j+1,k,kd)
                CALL gtvlr_bn(uh,dmngmtry,i,j,k+1,ke)
                CALL gt_distance(dmngmtry,i,j+1,k,'k-1',delta) !delta z
                k=k+1 !Artificio para utilizar subrutina get_distance de geometria
                CALL gt_distance(dmngmtry,i,j+1,k,'k-1',delte) !delta z
                k=k-1
                CALL gt_distance(dmngmtry,i,j,k+1,'j-1',delti) !delta y
                j=j+1 !Artificio para utilizar subrutina get_distance de geometria
                CALL gt_distance(dmngmtry,i,j,k+1,'j-1',delto) !delta y
                j=j-1
                IF ((acto/=0).AND.(acti/=0)) THEN 
                    Teq(a,b)=-(1/dy*(armncmn(kc,kd)/(delta+delte))+(1/dz*(armncmn(kc,ke)/(delti+delto))))
                ELSE IF (acto/=0) THEN
                    Teq(a,b)=-(1/dy*(armncmn(kc,kd)/(delta+delte)))
                ELSE IF (acti/=0) THEN
                    Teq(a,b)=-(1/dz*(armncmn(kc,ke)/(delti+delto)))
                ELSE
                    Teq(a,b)=0
                END IF     
            !Para bloque i,j+1,k-1   #12
            ELSE IF ((i==l).AND.(j==m-1).AND.(k==o+1)) THEN
                CALL gtvlr_bn(ua,dmngmtry,i,j+1,k+1,acto) 
                CALL gtvlr_bn(ua,dmngmtry,i,j-1,k-1,acti) 
                CALL gtvlr_bn(uh,dmngmtry,i,j+1,k,kd)
                CALL gtvlr_bn(uh,dmngmtry,i,j,k-1,ke)
                CALL gt_distance(dmngmtry,i,j+1,k,'k-1',delta) !delta z
                k=k+1 !Artificio para utilizar subrutina get_distance de geometria
                CALL gt_distance(dmngmtry,i,j+1,k,'k-1',delte) !delta z
                k=k-1
                CALL gt_distance(dmngmtry,i,j,k-1,'j-1',delti) !delta y
                j=j+1 !Artificio para utilizar subrutina get_distance de geometria
                CALL gt_distance(dmngmtry,i,j,k-1,'j-1',delto) !delta y
                j=j-1
                IF ((acto/=0).AND.(acti/=0)) THEN 
                    Teq(a,b)=(1/dy*(armncmn(kc,kd)/(delta+delte))+(1/dz*(armncmn(kc,ke)/(delti+delto))))
                ELSE IF (acto/=0) THEN
                    Teq(a,b)=(1/dy*(armncmn(kc,kd)/(delta+delte)))
                ELSE IF (acti/=0) THEN
                    Teq(a,b)=(1/dz*(armncmn(kc,ke)/(delti+delto)))
                ELSE
                    Teq(a,b)=0
                END IF  
            ELSE 
            END IF
        END DO anisocol
    END DO anisorow
    
    !
    CLOSE(u); CLOSE(ua); CLOSE(v); CLOSE(uhe)
    CLOSE(w); CLOSE(ul); CLOSE(ub)
    !
    WRITE (*,*) '*---|Finaliza construccion matriz densa con tensor anisotropo de conductividades|---*'
    !
END SUBROUTINE anisoTdns_bn