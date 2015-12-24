!!---------------------------------------------------------------------------------------------------------
!!Copyright (c) 2012, �lvarez-Villa O.D. & Universidad Polit�cnica de Valencia, All Right Reserved.
!!
!!This Program is distributed in the hope that they will be used but WITHOUT ANY WARRANTY. No autor or
!!distributor accepts responsability to anyone for consequences of using them or for whether they serve
!!any particula porpouse or work at all, unless he says so in writing. everyone is granted permission to
!!copy, modify an distribute this program, but only under the conditions that his notice and above
!!copyright notice remains intact. The author also requires citation in the use of the program for any
!!implementation as: "�lvarez-Villa, O.D (2012) SIMULACION EFICIENTE DE LAS RELACIONES RIO-ACUIFERO EN SISTEMAS 
!!DE UTILIZACION CONJUNTA MEDIANTE TECNICAS DE REDUCCION DE MODELOS LINEALES INVARIANTES EN EL TIEMPO.
!!Tesis Doctoral. Universidad Polit�cnica de Valencia".
!!---------------------------------------------------------------------------------------------------------
!!---------------------------------------------------------------------------------------------------------
!!Este m�dulo contiene subrutinas para simular el flujo en el acu�fero. Se obtienen los niveles
!!piezom�tricos en todos los bloques activos y los resultados se escriben en archivos de salida.
!!---------------------------------------------------------------------------------------------------------
!!
!MODULE fdsimulation
!    INTERFACE
!        SUBROUTINE simulation(ts,dt,prjctrt,prjctnm,tplgy,gmtry,steq,SF,topheconect,ijheconect, & 
!                              & qhcte,qhext,zns,isolv,itol,itmax,tol,prcndtntype,lfil,droptol,alp)
!            USE anisotype
!            IMPLICIT NONE
!            INTEGER(I4B), INTENT(INOUT)                          :: ts
!            REAL(DP), INTENT(IN)                                 :: dt
!            CHARACTER(*), INTENT(IN)                             :: prjctrt,prjctnm
!            TYPE(cnfnd_topology), INTENT(IN)                     :: tplgy
!            TYPE(confined_dp), INTENT(IN)                        :: gmtry
!            TYPE(sprs2_dp), INTENT(IN)                           :: steq
!            REAL(DP), DIMENSION(:), INTENT(INOUT), POINTER       :: SF 
!            INTEGER(I4B), DIMENSION(:), INTENT(INOUT), POINTER   :: topheconect
!            INTEGER(I4B), DIMENSION(:,:), INTENT(INOUT), POINTER :: ijheconect
!            REAL(DP), DIMENSION(:,:), INTENT(INOUT), POINTER     :: qhcte
!            REAL(DP), DIMENSION(:), INTENT(INOUT), POINTER       :: qhext   
!            INTEGER(I4B), DIMENSION(:), INTENT(INOUT), POINTER   :: zns
!            INTEGER(I4B), INTENT(IN)                             :: isolv,itol,itmax
!            REAL(DP), INTENT(IN)                                 :: tol    
!            INTEGER(I4B), INTENT(IN)                             :: prcndtntype,lfil
!            REAL(DP), INTENT(IN)                                 :: droptol,alp
!        END SUBROUTINE simulation
!    END INTERFACE
!    !!
!    INTERFACE
!        SUBROUTINE simulationstdy(prjctrt,prjctnm,tplgy,gmtry,steq,topheconect,ijheconect,qhcte, & 
!                                    & qhext,isolv,itol,itmax,tol,prcndtntype,lfil,droptol,alp)
!            USE anisotype
!            IMPLICIT NONE            !    
!            CHARACTER(*), INTENT(IN)                          :: prjctrt,prjctnm
!            TYPE(cnfnd_topology), INTENT(IN)                  :: tplgy
!            TYPE(confined_dp), INTENT(IN)                     :: gmtry
!            TYPE(sprs2_dp), INTENT(INOUT)                     :: steq
!            INTEGER(I4B), DIMENSION(:), INTENT(IN), POINTER   :: topheconect
!            INTEGER(I4B), DIMENSION(:,:), INTENT(IN), POINTER :: ijheconect
!            REAL(DP), DIMENSION(:,:), INTENT(IN), POINTER     :: qhcte
!            REAL(DP), DIMENSION(:), INTENT(IN), POINTER       :: qhext
!            INTEGER(I4B), INTENT(IN)                          :: isolv,itol,itmax
!            REAL(DP), INTENT(IN)                              :: tol    
!            INTEGER(I4B), INTENT(IN)                          :: prcndtntype,lfil
!            REAL(DP), INTENT(IN)                              :: droptol,alp
!        END SUBROUTINE simulationstdy
!    END INTERFACE
!    !!
!    INTERFACE
!        SUBROUTINE bldsystem(gmtry,tplgy,SF,topheconect,ijheconect,qhcte,qhext,htant,&
!                                  & dt,ar,ap,ahc,ahe,rowq,rowr,rowhe,rowhc,sa,b)
!            USE anisotype
!            IMPLICIT NONE
!            TYPE(confined_dp), INTENT(IN)                     :: gmtry
!            TYPE(cnfnd_topology), INTENT(IN)                  :: tplgy
!            REAL(DP), DIMENSION(:), INTENT(IN)                :: SF
!            INTEGER(I4B), DIMENSION(:), INTENT(IN), POINTER   :: topheconect
!            INTEGER(I4B), DIMENSION(:,:), INTENT(IN), POINTER :: ijheconect
!            REAL(DP), DIMENSION(:,:), INTENT(IN), POINTER     :: qhcte
!            REAL(DP), DIMENSION(:), INTENT(IN)                :: qhext
!            REAL(DP), DIMENSION(:), INTENT(IN)                :: htant
!            REAL(DP), INTENT(IN)                              :: dt
!            INTEGER(I4B), INTENT(IN)                          :: ar,ap,ahc,ahe
!            REAL(DP), DIMENSION(:) ,INTENT(IN)                :: rowq,rowr,rowhe,rowhc
!            TYPE(sprs2_dp), INTENT(INOUT)                     :: sa
!            REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT)  :: b
!        END SUBROUTINE bldsystem
!    END INTERFACE
!    !!
!    INTERFACE
!        SUBROUTINE bldsystemstdy(gmtry,tplgy,topheconect,ijheconect,qhcte,qhext, &
!                                  & ar,ap,ahc,ahe,rowq,rowr,rowhe,rowhc,sa,b)
!            USE anisotype
!            USE io,        ONLY : gtvlr
!            USE eigenutil, ONLY : eigenerror,assert_eq
!            !
!            IMPLICIT NONE
!            !
!            TYPE(confined_dp), INTENT(IN)                     :: gmtry
!            TYPE(cnfnd_topology), INTENT(IN)                  :: tplgy   
!            INTEGER(I4B), DIMENSION(:), INTENT(IN), POINTER   :: topheconect
!            INTEGER(I4B), DIMENSION(:,:), INTENT(IN), POINTER :: ijheconect
!            REAL(DP), DIMENSION(:,:), INTENT(IN), POINTER     :: qhcte
!            REAL(DP), DIMENSION(:), INTENT(IN)                :: qhext    
!            INTEGER(I4B), INTENT(IN)                          :: ar,ap,ahc,ahe
!            REAL(DP), DIMENSION(:) ,INTENT(IN)                :: rowq,rowr,rowhe,rowhc
!            TYPE(sprs2_dp), INTENT(IN)                        :: sa
!            REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT)  :: b
!        END SUBROUTINE bldsystemstdy
!    END INTERFACE
!    !!
!END MODULE
!!
!!***********************************************************************************************************************
!!******************RUTINA DE SIMULACI�N TRANSITORIA*********************************************************************
!!*********************@author: �lvarez-Villa O.D.***********************************************************************
!!***********************************************************************************************************************
!!
!!Subrutina que ejecuta la simulaci�n transitoria del flujo en el acu�fero.
!!Los param�tros de entrada son:
!!prjctrt        :Ruta del proyecto de trabajo.
!!prjnm          :Nombre del proyecto de trabajo.
!!tplgy          :La topologia del acuifero lineal
!!gmtry          :La geometria del acuifero
!!Teq            :La matriz de transmisividades equivalentes
!!topheconect    :Arreglo de indicadores de celdas activas donde existe conexi�n.
!!ijhcinect      :Arreglo que indican cuales celdas con CC de himpuesta se conectan con las activas
!!qhcte          :Arreglo con las transmisividades equivalentes de las celdas activas con conexi�n.
!!qhext          :Flujos generados por las H ext sobre las celdas activas
!!zns            :Arreglo con el n�mero de zonas de recarga, de condiciones de contorno y de pozos
!!Parametros para PBCG:
!!itol           :1,2,3,4, dependiendo del tipo de convergencia que se chequea
!!tol            :Tolerancia de convergencia deseada
!!itmax          :N�mero permitido de iteraciones
!!prcndtntype    :Entero que indica el tipo de precondicionamiento asignado
!!               :1: Precondicionamiento diagonal
!!               :2: Descomposici�n de Cholesky incompleta pic0
!!               :3: Descomposici�n LU incompleta sin relleno ilu0
!!               :4: Descomposici�n LU incompleta modificada sin relleno milu0
!!               :5: Descomposici�n LU incompleta con k nivel de relleno iluk
!!               :6: Descomposici�n LU incompleta con descarte simple ilud
!!               :7: Descomposici�n LU incompleta con estrategia doble de relleno
!!               :8: Precondicionador polinomial
!!lfil           :Nivel de llenado del precondicionamiento
!!droptol        :Valor para el descarte por tolerancia
!!alpha          :Nivel de conservaci�n por fila del precondicionamiento
!!
!SUBROUTINE simulation(ts,dt,prjctrt,prjctnm,tplgy,gmtry,steq,SF,topheconect,ijheconect, & 
!                        & qhcte,qhext,zns,isolv,itol,itmax,tol,prcndtntype,lfil,droptol,alp)
!    USE anisotype    
!    USE fdsimulation,     ONLY : bldsystem    
!    USE sprclinsolver,    ONLY : linpbcg, linpcg
!    USE eigenutil,        ONLY : eigenerror,assert_eq
!    USE io,               ONLY : gtvlr,openintensity,rdintensity,wrtht,prntmtrx
!    USE postprocess,      ONLY : rdheadswrt,rdqheswrt,rdvswrt,rdqhcswrt
!    USE preconditions,    ONLY : setprecond,cleanprecon
!    !
!    IMPLICIT NONE
!    !
!    INTEGER(I4B), INTENT(INOUT)                          :: ts
!    REAL(DP), INTENT(IN)                                 :: dt
!    CHARACTER(*), INTENT(IN)                             :: prjctrt,prjctnm
!    TYPE(cnfnd_topology), INTENT(IN)                     :: tplgy
!    TYPE(confined_dp), INTENT(IN)                        :: gmtry
!    TYPE(sprs2_dp), INTENT(INOUT)                        :: steq
!    REAL(DP), DIMENSION(:), INTENT(INOUT), POINTER       :: SF 
!    INTEGER(I4B), DIMENSION(:), INTENT(INOUT), POINTER   :: topheconect
!    INTEGER(I4B), DIMENSION(:,:), INTENT(INOUT), POINTER :: ijheconect
!    REAL(DP), DIMENSION(:,:), INTENT(INOUT), POINTER     :: qhcte
!    REAL(DP), DIMENSION(:), INTENT(INOUT), POINTER       :: qhext   
!    INTEGER(I4B), DIMENSION(:), INTENT(INOUT), POINTER   :: zns
!    INTEGER(I4B), INTENT(IN)                             :: isolv,itol,itmax
!    REAL(DP), INTENT(IN)                                 :: tol    
!    INTEGER(I4B), INTENT(IN)                             :: prcndtntype,lfil
!    REAL(DP), INTENT(IN)                                 :: droptol,alp            
!    !
!    TYPE(sprs2_dp)                      :: sa
!    TYPE(prcondition_dp)                :: prcndtn    
!    REAL(DP), DIMENSION(:), ALLOCATABLE :: rowq,rowr,rowhc,rowhe,ho,b,h
!    INTEGER(I4B)                        :: ur,uq,uhc,uhe,uri,uqi,uhci,uhei,uho,n
!    INTEGER(I4B)                        :: rzns,nwells,hczns,hezns
!    INTEGER(I4B)                        :: t1,t2,t3,t4,t,i,iter
!    CHARACTER(200)                      :: rtq,rthe,rthc,rtr,rtho
!    REAL(DP)                            :: val,err
!    !
!    PARAMETER(uho=56,uq=25,ur=93,uhc=64,uhe=89)
!    PARAMETER(uqi=26,uri=94,uhci=75,uhei=90)
!    !
!    !Verificar tama�o de los arreglos
!    n=assert_eq((/steq%n,SIZE(SF),SIZE(topheconect),tplgy%act/),'simulation')
!    !
!    !Abrir archivos de lectura directa
!    rtq=prjctrt//prjctnm//'.eigen.wel'                       !Archivo de localizaci�n de pozos
!    rtr=prjctrt//prjctnm//'.eigen.rcg'                       !Archivo de localizaci�n de zonas de recarga
!    rthc=prjctrt//prjctnm//'.eigen.qct'                      !Archivo de localizaci�n de CC de altura impuesta
!    rthe=prjctrt//prjctnm//'.eigen.qex'                      !Archivo de localizaci�n de CC de altura externa
!    rtho=prjctrt//prjctnm//'.eigen.inh'                      !Archivo de condiciones iniciales
!    OPEN(uq, FILE=rtq, STATUS='OLD', ACCESS='DIRECT', RECL=6, FORM='FORMATTED')
!    OPEN(ur, FILE=rtr, STATUS='OLD', ACCESS='DIRECT', RECL=6, FORM='FORMATTED')
!    OPEN(uhc, FILE=rthc, STATUS='OLD', ACCESS='DIRECT', RECL=6, FORM='FORMATTED')
!    OPEN(uhe, FILE=rthe, STATUS='OLD', ACCESS='DIRECT', RECL=6, FORM='FORMATTED')
!    OPEN(uho, FILE=rtho, STATUS='OLD', ACCESS='DIRECT', RECL=16, FORM='FORMATTED')
!    !
!    !Abrir archivos de lectura secuencial
!    CALL openintensity(prjctrt,prjctnm,'inr',uri,t1,rzns)     !Zonas homog�neas de recarga
!    CALL openintensity(prjctrt,prjctnm,'pmp',uqi,t2,nwells)   !Localizaci�n de pozos
!    CALL openintensity(prjctrt,prjctnm,'iqc',uhci,t3,hczns)   !Zonas de nivel impuesto
!    CALL openintensity(prjctrt,prjctnm,'iqe',uhei,t4,hezns)   !zonas de nivel externo
!    !
!    !Verificar zonas de archivo de ejecuci�n y de archivos de datos    
!    ts=assert_eq((/ts,t1,t2,t3,t4/),'simulation')             !Verificar horizonte de simulaci�n    
!    nwells=assert_eq((/nwells,zns(1)/),'simulation')          !Verificar n�mero de zonas
!    rzns=assert_eq((/rzns,zns(2)/),'simulation')
!    hezns=assert_eq((/hezns,zns(3)/),'simulation')
!    hczns=assert_eq((/hczns,zns(4)/),'simulation')
!    !
!    ALLOCATE(rowq(nwells),rowr(rzns),rowhc(hczns),rowhe(hezns),ho(n),b(n),h(n))    
!    !
!    ho=0.0_dp                                                 !Leer las alturas iniciales
!    DO i=1,n
!        CALL gtvlr(uho,gmtry%clmns,gmtry%rws,gmtry%lvls,tplgy%ijactv(i,1),tplgy%ijactv(i,2),1,ho(i))
!    END DO
!    t=0; h=ho   
!    CALL rdintensity(uhci,rowhc)                             !Lee las condiciones de contorno iniciales
!    CALL rdintensity(uhei,rowhe)                             !Escribe las alturas iniciales en un archivo
!    CALL wrtht(prjctrt,prjctnm,gmtry,tplgy,t,ho,rowhc,uhc)
!    !   
!    !Para el primer tiempo de simulaci�n se asume que h est� cerca de ho
!    !por lo cual al principio de las iteraciones del solver:
!    sa%n=steq%n; sa%len=steq%len
!    ALLOCATE(sa%val(steq%len),sa%irow(steq%len),sa%jcol(steq%len))
!    !    
!    PRINT*,'Inicia loop principal'
!    sim: DO t=1,ts                                           !Bucle de simulaci�n
!        PRINT*,'Tiempo',t
!        !        
!        CALL rdintensity(uqi,rowq)                           !Lee las filas de acciones externas 
!        CALL rdintensity(uri,rowr)                           !y condiciones de contorno
!        CALL rdintensity(uhci,rowhc)
!        CALL rdintensity(uhei,rowhe)
!        !       
!        PRINT*,'Construye sistema'                           !Llama el ensamblador del sistema
!        sa%val(1:sa%len)=steq%val(1:sa%len)
!        sa%irow(1:sa%len)=steq%irow(1:sa%len)
!        sa%jcol(1:sa%len)=steq%jcol(1:sa%len)        
!        iter=0; err=0.0_dp
!        CALL bldsystem(gmtry,tplgy,SF,topheconect,ijheconect,qhcte,qhext, & 
!                        & ho,dt,ur,uq,uhc,uhe,rowq,rowr,rowhe,rowhc,sa,b)
!        !        
!        PRINT*,'Resuelve sistema'                           !Construido el sistema, se soluciona por PBCG
!        CALL setprecond(sa,prcndtntype,lfil,droptol,alp,prcndtn)
!        SELECT CASE(isolv)
!            CASE(1)                
!                CALL linpcg(sa,prcndtn,b,h,tol,itmax,iter,err)
!            CASE(2)
!                IF(prcndtn%prcndtntype==2.OR.prcndtn%prcndtntype==8)&
!                            CALL eigenerror('unsoported precondining in: slvstdyh')
!                CALL linpbcg(sa,prcndtn,b,h,itol,tol,itmax,iter,err)
!            CASE DEFAULT
!                CALL eigenerror('illegal isolv in simulation in: slvstdyh')
!        END SELECT
!        CALL cleanprecon(prcndtn)
!        !        
!        PRINT*,'Escribe niveles piezom�tricos'              !Escribir alturas piezom�tricas simuladas
!        CALL wrtht(prjctrt,prjctnm,gmtry,tplgy,t,h,rowhc,uhc)
!        ho(1:n)=h(1:n)                                      !Iniciar alturas del intervalo siguiente
!        !
!    END DO sim
!    PRINT*,'Termin� simulaci�n'
!    !
!    !Destruir los arreglos con informaci�n para ahorrar espacio al postprocesamiento
!    DEALLOCATE(h,b,SF,rowq,rowr,rowhc,rowhe)                !Destruir las variables tipo    
!    DEALLOCATE(sa%val,sa%irow,sa%jcol,steq%val,steq%irow,steq%jcol)
!    !    
!    CLOSE(uq); CLOSE(ur); CLOSE(uhc)                        !Cerrar archivos de lectura
!    CLOSE(uhe); CLOSE(uho); CLOSE(uqi)
!    CLOSE(uri); CLOSE(uhci); CLOSE(uhei)
!    !
!    !Ahora se aplica el postprocesamiento para extraer las 
!    !variables de control de inter�s.
!    !
!    !1. Las alturas piezom�tricas
!    PRINT*,'Escritura de niveles piezom�tricos seleccionados'
!    CALL rdheadswrt(prjctrt,prjctnm,ts,gmtry)
!    !
!    !2. El volumen almacenado en zonas seleccionadas
!    PRINT*,'Escritura de volumenes almacenados en zonas seleccionadas'
!    CALL rdvswrt(prjctrt,prjctnm,ts,gmtry,tplgy)
!    !
!    !3. Caudal de intercambio con zonas de altura externa
!    PRINT*,'Escritura de caudales con r�os seleccionados'
!    CALL rdqheswrt(prjctrt,prjctnm,ts,gmtry,tplgy)
!    !
!    !4. Caudal de intercambio en zonas con nivel impuesto
!    PRINT*,'Escritura de caudales con niveles impuestos'
!    CALL rdqhcswrt(prjctrt,prjctnm,ts,gmtry,tplgy)
!    !
!END SUBROUTINE simulation
!!
!!***********************************************************************************************************************
!!******************RUTINA DE SIMULACI�N PERMANENTE**********************************************************************
!!*********************@author: �lvarez-Villa O.D.***********************************************************************
!!***********************************************************************************************************************
!!
!!prjctrt        :Ruta del proyecto de trabajo.
!!prjnm          :Nombre del proyecto de trabajo.
!!tplgy          :La topologia del acuifero lineal
!!gmtry          :La geometria del acuifero
!!Teq            :La matriz de transmisividades equivalentes
!!topheconect    :Arreglo de indicadores de celdas activas donde existe conexi�n.
!!ijhcinect      :Arreglo que indican cuales celdas con CC de himpuesta se conectan con las activas
!!qhcte          :Arreglo con las transmisividades equivalentes de las celdas activas con conexi�n.
!!qhext          :Flujos generados por las H ext sobre las celdas activas
!!zns            :Arreglo con el n�mero de zonas de recarga, de condiciones de contorno y de pozos
!!Parametros para PBCG:
!!itol           :1,2,3,4, dependiendo del tipo de convergencia que se chequea
!!tol            :Tolerancia de convergencia deseada
!!itmax          :N�mero permitido de iteraciones
!!prcndtntype    :Entero que indica el tipo de precondicionamiento asignado
!!               :1: Precondicionamiento diagonal
!!               :2: Descomposici�n de Cholesky incompleta pic0
!!               :3: Descomposici�n LU incompleta sin relleno ilu0
!!               :4: Descomposici�n LU incompleta modificada sin relleno milu0
!!               :5: Descomposici�n LU incompleta con k nivel de relleno iluk
!!               :6: Descomposici�n LU incompleta con descarte simple ilud
!!               :7: Descomposici�n LU incompleta con estrategia doble de relleno
!!               :8: Precondicionador polinomial
!!lfil           :Nivel de llenado del precondicionamiento
!!droptol        :Valor para el descarte por tolerancia
!!alpha          :Nivel de conservaci�n por fila del precondicionamiento
!!
!SUBROUTINE simulationstdy(prjctrt,prjctnm,tplgy,gmtry,steq,topheconect,ijheconect,qhcte, & 
!                          & qhext,isolv,itol,itmax,tol,prcndtntype,lfil,droptol,alp)
!    USE anisotype    
!    USE fdsimulation,     ONLY : bldsystemstdy    
!    USE sprclinsolver,    ONLY : linpbcg, linpcg
!    USE eigenutil,        ONLY : eigenerror,assert_eq
!    USE io,               ONLY : gtvlr,rdintensitysdty,wrthttsdy
!    USE stdypostprocess,  ONLY : rdheadswrtstdy,rdqheswrtstdy,rdvswrtstdy,rdqhcswrtstdy
!    USE preconditions,    ONLY : setprecond,cleanprecon
!    !
!    IMPLICIT NONE
!    !    
!    CHARACTER(*), INTENT(IN)                          :: prjctrt,prjctnm
!    TYPE(cnfnd_topology), INTENT(IN)                  :: tplgy
!    TYPE(confined_dp), INTENT(IN)                     :: gmtry
!    TYPE(sprs2_dp), INTENT(INOUT)                     :: steq
!    INTEGER(I4B), DIMENSION(:), INTENT(IN), POINTER   :: topheconect
!    INTEGER(I4B), DIMENSION(:,:), INTENT(IN), POINTER :: ijheconect
!    REAL(DP), DIMENSION(:,:), INTENT(IN), POINTER     :: qhcte
!    REAL(DP), DIMENSION(:), INTENT(IN), POINTER       :: qhext
!    INTEGER(I4B), INTENT(IN)                          :: isolv,itol,itmax
!    REAL(DP), INTENT(IN)                              :: tol    
!    INTEGER(I4B), INTENT(IN)                          :: prcndtntype,lfil
!    REAL(DP), INTENT(IN)                              :: droptol,alp
!    !
!    TYPE(sprs2_dp)                      :: sa
!    TYPE(prcondition_dp)                :: prcndtn    
!    REAL(DP), DIMENSION(:), ALLOCATABLE :: rowq,rowr,rowhc,rowhe,b,h
!    INTEGER(I4B)                        :: uexit,ur,uq,uhc,uhe,n,i,iter
!    INTEGER(I4B)                        :: rzns,nwells,hczns,hezns    
!    CHARACTER(200)                      :: rtq,rthe,rthc,rtr,rtexit
!    REAL(DP)                            :: val,err
!    !
!    PARAMETER(uq=25,ur=93,uhc=64,uhe=89,uexit=691)    
!    !
!    !Verificar tama�o de los arreglos
!    n=assert_eq((/steq%n,SIZE(topheconect),tplgy%act/),'simulationstdy')
!    !
!    !Abrir archivos de lectura directa
!    rtq=prjctrt//prjctnm//'.steady.wel'                      !Archivo de localizaci�n de pozos
!    rtr=prjctrt//prjctnm//'.steady.rcg'                      !Archivo de localizaci�n de zonas de recarga
!    rthc=prjctrt//prjctnm//'.steady.qct'                     !Archivo de localizaci�n de CC de altura impuesta
!    rthe=prjctrt//prjctnm//'.steady.qex'                     !Archivo de localizaci�n de CC de altura externa
!    rtexit=prjctrt//'Results\'//prjctnm//'.steady.out'       !Archivo de saludas estacionarias    
!    OPEN(uq, FILE=rtq, STATUS='OLD', ACCESS='DIRECT', RECL=6, FORM='FORMATTED')
!    OPEN(ur, FILE=rtr, STATUS='OLD', ACCESS='DIRECT', RECL=6, FORM='FORMATTED')
!    OPEN(uhc, FILE=rthc, STATUS='OLD', ACCESS='DIRECT', RECL=6, FORM='FORMATTED')
!    OPEN(uhe, FILE=rthe, STATUS='OLD', ACCESS='DIRECT', RECL=6, FORM='FORMATTED')
!    !
!    ALLOCATE(b(n),h(n))
!    !
!    !Abrir archivos de lectura secuencial
!    CALL rdintensitysdty(prjctrt,prjctnm,rowq,rowr,rowhc,rowhe)!Intensidades de las AE    
!    !
!    !Verificar zonas de archivo de ejecuci�n y de archivos de datos    
!    nwells=SIZE(rowq); rzns=SIZE(rowr)
!    hezns=SIZE(rowhe); hczns=SIZE(rowhc)
!    !    
!    PRINT*,'Construye sistema'                                 !Llama el ensamblador del sistema    
!    sa%n=steq%n; sa%len=steq%len
!    ALLOCATE(sa%val(steq%len),sa%irow(steq%len),sa%jcol(steq%len))
!    sa%val(1:sa%len)=steq%val(1:sa%len)
!    sa%irow(1:sa%len)=steq%irow(1:sa%len)
!    sa%jcol(1:sa%len)=steq%jcol(1:sa%len)        
!    iter=0; err=0.0_dp; h=0.0_dp
!    CALL bldsystemstdy(gmtry,tplgy,topheconect,ijheconect,qhcte,qhext, & 
!                       & ur,uq,uhc,uhe,rowq,rowr,rowhe,rowhc,sa,b)
!    !        
!    PRINT*,'Resuelve sistema'                                 !Construido el sistema, se soluciona por PBCG
!!    CALL setprecond(sa,prcndtntype,lfil,droptol,alp,prcndtn)
!!    SELECT CASE(isolv)
!!        CASE(1)                
!!            CALL linpcg(sa,prcndtn,b,h,tol,itmax,iter,err)
!!        CASE(2)
!!            IF(prcndtn%prcndtntype==2.OR.prcndtn%prcndtntype==8)&
!!                        CALL eigenerror('unsoported precondining in: slvstdyh')
!!            CALL linpbcg(sa,prcndtn,b,h,itol,tol,itmax,iter,err)
!!        CASE DEFAULT
!!            CALL eigenerror('illegal isolv in simulation in: slvstdyh')
!!    END SELECT
!!    CALL cleanprecon(prcndtn)
!    !        
!    PRINT*,'Escribe niveles piezom�tricos'                    !Escribir alturas piezom�tricas simuladas
!    CALL wrthttsdy(prjctrt,prjctnm,gmtry,tplgy,h,rowhc,uhc)
!    OPEN(uexit,FILE=rtexit,STATUS='REPLACE')
!    !
!!    !Destruir los arreglos con informaci�n para ahorrar espacio al postprocesamiento
!!    DEALLOCATE(h,b)                                           !Destruir las variables tipo    
!!    DEALLOCATE(sa%val,sa%irow,sa%jcol,steq%val,steq%irow,steq%jcol)
!!    !
!!    !Ahora se aplica el postprocesamiento para extraer las 
!!    !variables de control de inter�s.
!!    !
!!    !1. Las alturas piezom�tricas
!!    PRINT*,'Escritura de niveles piezom�tricos seleccionados'
!!    CALL rdheadswrtstdy(uexit,prjctrt,prjctnm,gmtry)
!!    !
!!    !2. El volumen almacenado en zonas seleccionadas
!!    PRINT*,'Escritura de volumenes almacenados en zonas seleccionadas'
!!    CALL rdvswrtstdy(uexit,prjctrt,prjctnm,gmtry,tplgy)
!!    !
!!    !3. Caudal de intercambio con zonas de altura externa
!!    PRINT*,'Escritura de caudales con r�os seleccionados'
!!    CALL rdqheswrtstdy(uexit,prjctrt,prjctnm,gmtry,tplgy,rowhe)
!!    !
!!    !4. Caudal de intercambio en zonas con nivel impuesto
!!    PRINT*,'Escritura de caudales con niveles impuestos'
!!    CALL rdqhcswrtstdy(uexit,prjctrt,prjctnm,gmtry,tplgy)
!!    !
!!    DEALLOCATE(rowq,rowr,rowhc,rowhe)
!!    CLOSE(uq); CLOSE(ur)                                      !Cerrar archivos de lectura
!!    CLOSE(uhe); CLOSE(uhc); CLOSE(uexit)
!    !
!END SUBROUTINE simulationstdy
!!
!!***********************************************************************************************************************
!!***************RUTINA DE ENSAMBLE DE SISTEMAS PARA DT******************************************************************
!!*********************@author: �lvarez-Villa O.D.***********************************************************************
!!***********************************************************************************************************************
!!
!!Subrutina que construye el sistema de acuaciones lineales a resolver
!!en cada intervalo temporal. Los par�metros son:
!!tplgy          :La topologia del acuifero lineal
!!gmtry          :La geometria del acuifero
!!Teq            :La matriz de transmisividades equivalentes
!!topheconect    :Arreglo de indicadores de celdas activas donde existe conexi�n.
!!ijhcinect      :Arreglo que indican cuales celdas con CC de himpuesta se conectan con las activas
!!qhcte          :Arreglo con las transmisividades equivalentes de las celdas activas con conexi�n.
!!qhext          :Flujos generados por las H ext sobre las celdas activas.
!!hant           :Alturas piezom�tricas obtenidas en el intervalo anterior.
!!dt             :Duraci�n en d�as del intervalo temporal.
!!ar,ap,ahc,ahe  :Localizaci�n en memoria de los archivos formateados que indican la existencia de AE.
!!rowq           :Bombeos por pozos, para el tiempo en consideraci�n
!!rowr           :Recargas por zonas, para el tiempo en consideraci�n
!!rowhe          :Alturas para CC de nivel externo, para el tiempo en consideraci�n
!!rowhc          :Alturas para CC de nivel impuesto, para el tiempo en consideraci�n
!!                de las acciones exteriores, r:recarga, p:bonbeo,hc:Niveles impuestos, he:Nivel exterior
!!Las salidas de la subrutina son:
!!sa             :Matriz dispersa de coeficientes.
!!b              :Vector de t�rminos independientes.
!!
!SUBROUTINE bldsystem(gmtry,tplgy,SF,topheconect,ijheconect,qhcte,qhext,htant,&
!                        & dt,ar,ap,ahc,ahe,rowq,rowr,rowhe,rowhc,sa,b)
!    USE anisotype
!    USE io,        ONLY : gtvlr
!    USE eigenutil, ONLY : eigenerror,assert_eq
!    !
!    IMPLICIT NONE
!    !
!    TYPE(confined_dp), INTENT(IN)                     :: gmtry
!    TYPE(cnfnd_topology), INTENT(IN)                  :: tplgy
!    REAL(DP), DIMENSION(:), INTENT(IN)                :: SF
!    INTEGER(I4B), DIMENSION(:), INTENT(IN), POINTER   :: topheconect
!    INTEGER(I4B), DIMENSION(:,:), INTENT(IN), POINTER :: ijheconect
!    REAL(DP), DIMENSION(:,:), INTENT(IN), POINTER     :: qhcte
!    REAL(DP), DIMENSION(:), INTENT(IN)                :: qhext
!    REAL(DP), DIMENSION(:), INTENT(IN)                :: htant
!    REAL(DP), INTENT(IN)                              :: dt
!    INTEGER(I4B), INTENT(IN)                          :: ar,ap,ahc,ahe
!    REAL(DP), DIMENSION(:) ,INTENT(IN)                :: rowq,rowr,rowhe,rowhc
!    TYPE(sprs2_dp), INTENT(INOUT)                     :: sa
!    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT)  :: b
!    !
!    INTEGER(I4B)                        :: i,j,k,l,a,n
!    INTEGER(I4B)                        :: cont,znwell,znrcg,znhc,znhe
!    REAL(DP), DIMENSION(:), ALLOCATABLE :: q
!    !                                                  !Verificar tama�o de los arreglos    
!    n=assert_eq((/SIZE(SF),SIZE(htant),SIZE(topheconect),tplgy%act,sa%n/),'bldsystem')
!    ALLOCATE(b(n),q(n))                                !Formo sa: sparce([Teq]-[SF]/dt)
!    b=0.0_dp; q=0.0_dp; cont=0 
!    !  
!    itdiagsa: DO k=1,sa%len                            !Reemplazar en la diagonal de sa
!        !
!        IF( sa%irow(k).EQ.sa%jcol(k) )THEN 
!            cont=cont+1
!            b(cont)=(-1.0_dp)*SF(cont)/dt
!            sa%val(k)=sa%val(k)+b(cont)
!            b(cont)=b(cont)*htant(cont)
!            
!            !Ahora se modifican las entradas de q
!            i=tplgy%ijactv(cont,1); j=tplgy%ijactv(cont,2)            
!            !
!            !1. Los bombeos
!            !identificar si en esta celda activa hay pozo y su convenci�n
!            CALL gtvlr(ap,gmtry%clmns,gmtry%rws,gmtry%lvls,i,j,1,znwell)
!            IF(znwell.NE.0)THEN   !Por convenci�n: Bombeos negativos
!               q(cont)=q(cont)-rowq(znwell)   
!            END IF
!            !WRITE(*,*)'i',i,'j',j,'zona=',znwell,'recarga=',rowq(znwell)            
!            !
!            !2. Las recargas
!            !identificar cual es la zona de recarga de esta celda activa
!            CALL gtvlr(ar,gmtry%clmns,gmtry%rws,gmtry%lvls,i,j,1,znrcg)
!            IF(znrcg.NE.0)THEN   !Por convenci�n: Recargas positivas
!               q(cont)=q(cont)+rowr(znrcg)*( gmtry%X(j+1)-gmtry%X(j) )*( gmtry%Y(i+1)-gmtry%Y(i) )   !R=r*A
!            END IF
!            !WRITE(*,*)'i',i,'j',j,'zona=',znrcg,'recarga=',rowr(znrcg)
!            !
!            !3. Las condiciones de contorno de nivel externo
!            !identificar la zona a la que pertenece la celda activa
!            CALL gtvlr(ahe,gmtry%clmns,gmtry%rws,gmtry%lvls,i,j,1,znhe)
!            IF(znhe.NE.0)THEN   !Por convenci�n: Recargas positivas
!               q(cont)=q(cont)+rowhe(znhe)*qhext(cont)   !q=T*he
!            END IF
!            !
!            !4. Las condiciones de contorno de nivel impuesto
!            !identifica si la celda activa tiene conexi�n
!            IF(topheconect(cont).NE.0)THEN
!                !si tiene conexi�n identifica, realiza la acumulaci�n de caudales                
!                IF( ijheconect(topheconect(cont),1).EQ.1 )THEN   
!                    a=i-1; l=j                          !Bloque inferior
!                    CALL gtvlr(ahc,gmtry%clmns,gmtry%rws,gmtry%lvls,a,l,1,znhc)
!                    q(cont)=q(cont)+rowhc(znhc)*qhcte(topheconect(cont),1)
!                END IF
!                IF( ijheconect(topheconect(cont),2).EQ.1 )THEN   
!                    a=i+1; l=j                          !Bloque superior
!                    CALL gtvlr(ahc,gmtry%clmns,gmtry%rws,gmtry%lvls,a,l,1,znhc)
!                    q(cont)=q(cont)+rowhc(znhc)*qhcte(topheconect(cont),2)
!                END IF
!                IF( ijheconect(topheconect(cont),3).EQ.1 )THEN   
!                    a=i; l=j-1                          !Bloque izquierdo
!                    CALL gtvlr(ahc,gmtry%clmns,gmtry%rws,gmtry%lvls,a,l,1,znhc)
!                    q(cont)=q(cont)+rowhc(znhc)*qhcte(topheconect(cont),3)
!                END IF
!                IF( ijheconect(topheconect(cont),4).EQ.1 )THEN   
!                    a=i; l=j+1                          !Bloque derecho
!                    CALL gtvlr(ahc,gmtry%clmns,gmtry%rws,gmtry%lvls,a,l,1,znhc)
!                    q(cont)=q(cont)+rowhc(znhc)*qhcte(topheconect(cont),4)
!                END IF
!            END IF
!            !            
!            b(cont)=b(cont)-q(cont)                     !Vector de t�rminos independientes
!        END IF
!        !
!    END DO itdiagsa    
!    !
!END SUBROUTINE bldsystem
!!
!!***********************************************************************************************************************
!!***************RUTINA DE ENSAMBLE DE SISTEMAS PERMANENTES**************************************************************
!!*********************@author: �lvarez-Villa O.D.***********************************************************************
!!***********************************************************************************************************************
!!Subrutina que construye el sistema de acuaciones lineales a resolver
!!para un modelo permanente. Los par�metros son:
!!tplgy          :La topologia del acuifero lineal
!!gmtry          :La geometria del acuifero
!!Teq            :La matriz de transmisividades equivalentes
!!topheconect    :Arreglo de indicadores de celdas activas donde existe conexi�n.
!!ijhcinect      :Arreglo que indican cuales celdas con CC de himpuesta se conectan con las activas
!!qhcte          :Arreglo con las transmisividades equivalentes de las celdas activas con conexi�n.
!!qhext          :Flujos generados por las H ext sobre las celdas activas.
!!hant           :Alturas piezom�tricas obtenidas en el intervalo anterior.
!!dt             :Duraci�n en d�as del intervalo temporal.
!!ar,ap,ahc,ahe  :Localizaci�n en memoria de los archivos formateados que indican la existencia de AE.
!!rowq           :Bombeos por pozos, para el tiempo en consideraci�n
!!rowr           :Recargas por zonas, para el tiempo en consideraci�n
!!rowhe          :Alturas para CC de nivel externo, para el tiempo en consideraci�n
!!rowhc          :Alturas para CC de nivel impuesto, para el tiempo en consideraci�n
!!                de las acciones exteriores, r:recarga, p:bonbeo,hc:Niveles impuestos, he:Nivel exterior
!!Las salidas de la subrutina son:
!!sa             :Matriz dispersa de coeficientes.
!!b              :Vector de t�rminos independientes.
!!
!SUBROUTINE bldsystemstdy(gmtry,tplgy,topheconect,ijheconect,qhcte,qhext, &
!                          & ar,ap,ahc,ahe,rowq,rowr,rowhe,rowhc,sa,b)
!    USE anisotype
!    USE io,        ONLY : gtvlr
!    USE eigenutil, ONLY : eigenerror,assert_eq
!    !
!    IMPLICIT NONE
!    !
!    TYPE(confined_dp), INTENT(IN)                     :: gmtry
!    TYPE(cnfnd_topology), INTENT(IN)                  :: tplgy   
!    INTEGER(I4B), DIMENSION(:), INTENT(IN), POINTER   :: topheconect
!    INTEGER(I4B), DIMENSION(:,:), INTENT(IN), POINTER :: ijheconect
!    REAL(DP), DIMENSION(:,:), INTENT(IN), POINTER     :: qhcte
!    REAL(DP), DIMENSION(:), INTENT(IN)                :: qhext    
!    INTEGER(I4B), INTENT(IN)                          :: ar,ap,ahc,ahe
!    REAL(DP), DIMENSION(:) ,INTENT(IN)                :: rowq,rowr,rowhe,rowhc
!    TYPE(sprs2_dp), INTENT(IN)                        :: sa
!    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT)  :: b
!    !
!    REAL(DP)                            :: x
!    INTEGER(I4B)                        :: i,j,k,l,a,n
!    INTEGER(I4B)                        :: cont,znwell,znrcg,znhc,znhe
!    REAL(DP), DIMENSION(:), ALLOCATABLE :: q
!    !                                                  !Verificar tama�o de los arreglos    
!    n=assert_eq((/SIZE(topheconect),tplgy%act,sa%n/),'bldsystem')
!    ALLOCATE(b(n),q(n))                                !Formo sa: sparce([Teq]-[SF]/dt)
!    b=0.0_dp; q=0.0_dp; cont=0 
!    !  
!    itdiagsa: DO k=1,sa%len                            !Reemplazar en la diagonal de sa
!        !
!         IF( sa%irow(k).EQ.sa%jcol(k) )THEN 
!            cont=cont+1
!            !
!            !Ahora se modifican las entradas de q
!            i=tplgy%ijactv(cont,1); j=tplgy%ijactv(cont,2)
!            !
!            !1. Los bombeos
!            !identificar si en esta celda activa hay pozo y su convenci�n
!            CALL gtvlr(ap,gmtry%clmns,gmtry%rws,gmtry%lvls,i,j,1,znwell)
!            IF(znwell.NE.0)THEN   !Por convenci�n: Bombeos negativos
!                x=rowq(znwell) 
!               q(cont)=q(cont)-x
!            END IF
!            !
!            !2. Las recargas
!            !identificar cual es la zona de recarga de esta celda activa
!            CALL gtvlr(ar,gmtry%clmns,gmtry%rws,gmtry%lvls,i,j,1,znrcg)
!            IF(znrcg.NE.0)THEN   !Por convenci�n: Recargas positivas
!                x=rowr(znrcg)*( gmtry%X(j+1)-gmtry%X(j) )*( gmtry%Y(i+1)-gmtry%Y(i) )
!               q(cont)=q(cont)+ x  !R=r*A
!            END IF
!            !
!            !3. Las condiciones de contorno de nivel externo
!            !identificar la zona a la que pertenece la celda activa
!            CALL gtvlr(ahe,gmtry%clmns,gmtry%rws,gmtry%lvls,i,j,1,znhe)
!            IF(znhe.NE.0)THEN   !Por convenci�n: Recargas positivas
!               x=rowhe(znhe)*qhext(cont)
!               q(cont)=q(cont)+x  !q=T*he
!            END IF
!            !
!            !4. Las condiciones de contorno de nivel impuesto
!            !identifica si la celda activa tiene conexi�n
!            IF(topheconect(cont).NE.0)THEN
!                !si tiene conexi�n identifica, realiza la acumulaci�n de caudales                
!                IF( ijheconect(topheconect(cont),1).EQ.1 )THEN   
!                    a=i-1; l=j                          !Bloque inferior
!                    CALL gtvlr(ahc,gmtry%clmns,gmtry%rws,gmtry%lvls,a,l,1,znhc)
!                    q(cont)=q(cont)+rowhc(znhc)*qhcte(topheconect(cont),1)
!                END IF
!                IF( ijheconect(topheconect(cont),2).EQ.1 )THEN   
!                    a=i+1; l=j                          !Bloque superior
!                    CALL gtvlr(ahc,gmtry%clmns,gmtry%rws,gmtry%lvls,a,l,1,znhc)
!                    q(cont)=q(cont)+rowhc(znhc)*qhcte(topheconect(cont),2)
!                END IF
!                IF( ijheconect(topheconect(cont),3).EQ.1 )THEN   
!                    a=i; l=j-1                          !Bloque izquierdo
!                    CALL gtvlr(ahc,gmtry%clmns,gmtry%rws,gmtry%lvls,a,l,1,znhc)
!                    q(cont)=q(cont)+rowhc(znhc)*qhcte(topheconect(cont),3)
!                END IF
!                IF( ijheconect(topheconect(cont),4).EQ.1 )THEN   
!                    a=i; l=j+1                          !Bloque derecho
!                    CALL gtvlr(ahc,gmtry%clmns,gmtry%rws,gmtry%lvls,a,l,1,znhc)
!                    q(cont)=q(cont)+rowhc(znhc)*qhcte(topheconect(cont),4)
!                END IF
!            END IF
!            !            
!            b(cont)=b(cont)-q(cont)                     !Vector de t�rminos independientes
!            !
!         END IF
!         !
!    END DO itdiagsa
!    !
!END SUBROUTINE bldsystemstdy