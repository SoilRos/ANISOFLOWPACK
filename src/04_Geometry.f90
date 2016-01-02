!-----------------------------------------------------------------------------------------------------------------------------------
!Copyright (c) 2013, Alvarez-Villa O.D., Alvarez-Rendon J.P., Perez Kevin  GOTTA Ingenieria SAS, All Right Reserved.
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
MODULE geometry
    !!
    INTERFACE bld_domainff
        SUBROUTINE bld_domainff_dp(prjctrt,flnm,X,Y,Z,dmngmtry)    
            USE anisotype
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN)                                 :: prjctrt,flnm
            TYPE(domain), INTENT(OUT)                                :: dmngmtry
            REAL(DP), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(OUT) :: X,Y,Z
        END SUBROUTINE bld_domainff_dp
    END INTERFACE bld_domainff
    !!
    INTERFACE bld_domain
        SUBROUTINE bld_domain_dp(X,Y,Z,rws,clmns,lvls,dmngmtry)
            USE anisotype
            IMPLICIT NONE
            INTEGER(I4B), INTENT(IN)                                    :: rws,clmns,lvls
            REAL(DP), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(IN)     :: X,Y,Z
            TYPE(domain), INTENT(OUT)                                  :: dmngmtry
        END SUBROUTINE bld_domain_dp
    END INTERFACE bld_domain
    !!
    INTERFACE 
        SUBROUTINE dstry_domain(dmngmtry,switch)    
            USE anisotype
            IMPLICIT NONE
            LOGICAL(LGT), INTENT(IN)    :: switch
            TYPE(domain), INTENT(INOUT) :: dmngmtry
        END SUBROUTINE dstry_domain
    END INTERFACE
    !!7
    INTERFACE gt_vxlcntr
        SUBROUTINE gt_vxlcntr_dp(dmngmtry,i,j,k,x,y,z)
            USE anisotype
            IMPLICIT NONE
            INTEGER(I4B), INTENT(IN)    :: i,j,k
            REAL(DP), INTENT(OUT)       :: x,y,z
            TYPE(domain), INTENT(IN)   :: dmngmtry
        END SUBROUTINE gt_vxlcntr_dp
    END INTERFACE gt_vxlcntr
    !!
    INTERFACE gt_dmnbnds
        SUBROUTINE gt_dmnbnds_dp(dmngmtry,xi,yi,zi,xf,yf,zf)
            USE anisotype
            IMPLICIT NONE
            REAL(DP), INTENT(OUT)       :: xi,yi,zi,xf,yf,zf
            TYPE(domain), INTENT(IN)   :: dmngmtry
        END SUBROUTINE gt_dmnbnds_dp
    END INTERFACE gt_dmnbnds
    
    INTERFACE gt_delta
        SUBROUTINE gt_delta_dp(dmngmtry,i,j,k,dx,dy,dz)
            USE anisotype
            IMPLICIT NONE        
            INTEGER(I4B), INTENT(IN)    :: i,j,k
            REAL(DP), INTENT(OUT)       :: dx,dy,dz
            TYPE(domain), INTENT(IN)   :: dmngmtry
        END SUBROUTINE gt_delta_dp
    END INTERFACE gt_delta
    
    INTERFACE gt_distance
        SUBROUTINE gt_distance_dp(dmngmtry,i,j,k,c,dist)
            USE anisotype
            IMPLICIT NONE 
            INTEGER(I4B), INTENT(IN)    :: i,j,k
            CHARACTER(3), INTENT(IN)    :: c
            REAL(DP), INTENT(OUT)       :: dist
            TYPE(domain), INTENT(IN)    :: dmngmtry
        END SUBROUTINE gt_distance_dp
    END INTERFACE gt_distance
END MODULE geometry
!

!***********************************************************************************************************************************
!*************************** CONSTRUCCION DE UNA VARIABLE TIPO GEOMETRIA DESDE UN ARCHIVO EN DISCO *********************************
!*****************************@authores: Alvarez-Villa O.D. & Rendon-Alvarez J.P. & Perez K. ***************************************
!***********************************************************************************************************************************
!
!DESCRIPCION
!NOMBRE: bld_domainff (build domain from file).
!Subrutina que construye la variable tipo domain (geometria) a partir de una archivo de geometria almacenado en disco. Para tal fin
!utiliza la subrutina bld_domain, que asigna a los punteros caracteristicos de la variable tipo las variables alamacenadas en RAM
!a partir de la lectura realizada.

!DATOS DE ENTRADA:
!prjctrt        : Ruta del proyecto de trabajo                                                                            {caracter}
!flnm           : Nombre del archivo de texto en el disco duro con su extencion                                           {caracter}
!X,Y,Z          : Arreglos donde se almacenaron las coordenadas                                                         {real}(:)[O]
!
!DATOS DE SALIDA:
!dmngmtry       : Variable tipo de geometria                                                                                {domain}
!
!IMPORTANTE:
!i=Contador de columnas (x), j=Contador de filas (y), k=Contador de niveles (z) 

SUBROUTINE bld_domainff_dp(prjctrt,flnm,X,Y,Z,dmngmtry)    
    USE anisotype
    USE util,       ONLY : error
    USE geometry,   ONLY : bld_domain
    !
    IMPLICIT NONE
#include <petsc/finclude/petscsys.h>
    PetscErrorCode                                           :: ierr
    PetscMPIInt                                              :: tprocess, iprocess,status(MPI_STATUS_SIZE)
    !
    CHARACTER(*), INTENT(IN)                                 :: prjctrt,flnm
    TYPE(domain), INTENT(OUT)                                :: dmngmtry
    REAL(DP), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(OUT) :: X,Y,Z
    !
    CHARACTER(1)                                :: enter
    CHARACTER(2)                                :: rtdlmtng
    CHARACTER(8)                                :: header
    CHARACTER(200)                              :: rt        
    INTEGER(I4B)                                :: u,i,j,k,clmns,rws,lvls,hrs,mnts,a
    REAL(DP)                                    :: timea,timeb,exetime,scnds
    LOGICAL(LGT)                                :: flexist
    !
    PARAMETER(u=29)
    !
    CALL MPI_comm_rank(MPI_COMM_WORLD,iprocess,ierr)
    CALL MPI_comm_size(MPI_COMM_WORLD,tprocess,ierr)
    CALL MPI_Barrier(MPI_COMM_WORLD,ierr)                                                                                               !Barrier: Sincroniza todos los procesadores al inicio, esto para ubicar más fácil posibles errores.
    IF (iprocess==0) THEN
        WRITE(*,*) '*---|Importando archivo de geometria|---*'                                                                     
        !Definicion del caracter de delimitacion de rutas de acuerdo con el sistema operativo 
        IF (INDEX(prjctrt,'/')==0) THEN
            rtdlmtng='\\'                                                                                                               !Delimitador en MS-DOS
        ELSE
            rtdlmtng='//'                                                                                                               !Delimitador en Unix-OS
        END IF
        !Verificacion de la existencia y pertinencia de los datos de entrada
        IF (LEN_TRIM(TRIM(ADJUSTL(prjctrt))//rtdlmtng//TRIM(ADJUSTL(flnm))) > 200) THEN
            CALL error('En bld_domainff_dp: La ruta del archivo geometria es demasiado larga')    
        ELSE
            rt=TRIM(ADJUSTL(prjctrt))//rtdlmtng//TRIM(ADJUSTL(flnm))                                                                    !Archivo en el que  se almacena el dominio espacial        
        END IF
        INQUIRE(FILE=TRIM(ADJUSTL(rt)),EXIST=flexist)
        IF (.NOT.flexist) CALL error('En bld_domainff_dp: El archivo geometria no existe')                                              !Verifica que el archivo con el dominio exista y est� en la ruta correcta
        !Lectura del archivo que contiene el dominio espacial
        OPEN(u,FILE=rt,STATUS='OLD',ACTION='READ')                                                                                      !Abre el archivo del dominio espacial
        READ(u,'(A8,I10,A1)')header,clmns,enter                                                                                         !Numero de columnas
        READ(u,'(A8,I10,A1)')header,rws,enter                                                                                           !Numero de filas
        READ(u,'(A8,I10,A1)')header,lvls,enter                                                                                          !Numero de niveles
        ALLOCATE(X(clmns+1),Y(rws+1),Z(lvls+1))                                                                                         !Aloja los arreglos de coordenadas en precisi�n doble
        DO i=1,clmns+1
            READ(u,'(SP,ES18.11E2,A1)')X(i),enter                                                                                       !Coordenadas en la direccion X
        END DO
        DO j=1,rws+1
            READ(u,'(SP,ES18.11E2,A1)')Y(j),enter                                                                                       !Coordenadas en la direccion Y
        END DO
        DO k=1,lvls+1
            READ(u,'(SP,ES18.11E2,A1)')Z(k),enter                                                                                       !Coordenadas en la direccion Z
        END DO
        CLOSE(u)
        CALL bld_domain(X,Y,Z,rws,clmns,lvls,dmngmtry)
        
        ! Envia a los otros procesadores la cantidad de niveles, filas y columnas.
        DO a=1,tprocess-1
            CALL MPI_SEND(dmngmtry%lvls ,1,MPI_INTEGER,a,200,MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(dmngmtry%rws  ,1,MPI_INTEGER,a,201,MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(dmngmtry%clmns,1,MPI_INTEGER,a,202,MPI_COMM_WORLD,ierr)
        END DO
        CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
        WRITE(*,*) '*---|La creacion de la variable tipo geometria ha sido exitosa|---*'                                                !Imprime el mensaje de finalizacion satisfactoria del programa
    ELSE
        CALL MPI_RECV(dmngmtry%lvls ,1,MPI_INTEGER,MPI_ANY_SOURCE,200,MPI_COMM_WORLD, status, ierr)
        CALL MPI_RECV(dmngmtry%rws  ,1,MPI_INTEGER,MPI_ANY_SOURCE,201,MPI_COMM_WORLD, status, ierr)
        CALL MPI_RECV(dmngmtry%clmns,1,MPI_INTEGER,MPI_ANY_SOURCE,202,MPI_COMM_WORLD, status, ierr)
        CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
    END IF
    !
END SUBROUTINE bld_domainff_dp
!!
!***********************************************************************************************************************************
!************************** ASGINACION DE ATRIBUTOS ALAMACENADOS EN RAM A VARIABLE TIPO GEOMETRIA **********************************
!*****************************@authores: Alvarez-Villa O.D. & Rendon-Alvarez J.P. & Perez K. ***************************************
!***********************************************************************************************************************************
!
!DESCRIPCION:
!Nombre: bld_domain(build domain)
!Subrutina auxiliar de bld_domainff. Asigna a los atributos de la variable tipo geometria las variables cargadas en RAM por 
!bld_domainff.
!
!DATOS DE ENTRADA:
!X          : Arreglo de coordenadas en la direccion X                                                                  {real}(:)[O]
!Y          : Arreglo de coordenadas en la direccion Y                                                                  {real}(:)[O]
!Z          : Arreglo de coordenadas en la direccion Z                                                                  {real}(:)[O]
!clmns      : Numero de columnas                                                                                            {entero}
!rws        : Numero de filas                                                                                               {entero}
!lvls       : Numero de niveles                                                                                             {entero}
!
!DATOS DE SALIDA:
!dmngmtry   : Variable tipo dominio espacial                                                                                {domain}
SUBROUTINE bld_domain_dp(X,Y,Z,rws,clmns,lvls,dmngmtry)
    USE anisotype
    USE util, ONLY : error, assert_eq
    !
    IMPLICIT NONE
    !
    INTEGER(I4B), INTENT(IN)                                    :: rws,clmns,lvls
    REAL(DP), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(IN)     :: X,Y,Z
    TYPE(domain), INTENT(OUT)                                   :: dmngmtry
    !
    !Verificaci�n de consistencia en los datos de entrada
    IF (SIZE(X,1)/=clmns+1) CALL error('En bld_domain_sp: Dimensiones incongruentes del vector x')
    IF (SIZE(Y,1)/=rws+1)   CALL error('En bld_domain_sp: Dimensiones incongruentes del vector y')
    IF (SIZE(Z,1)/=lvls+1)  CALL error('En bld_domain_sp: Dimensiones incongruentes del vector z')
    !Construcci�n del dominio espacial
    dmngmtry%exist=.true.                                                                                                           !Garantiza la existencia de la variable tipo
    dmngmtry%infile=.false.                                                                                                         !no daTa ya que la variable esta en RAM
    dmngmtry%flind=-9999                                                                                                            !...todos los par�metros de almacenamiento en...
    dmngmtry%rtdom=''                                                                                                               !Atributo reservado para hacer la lEctura directa de disco
    dmngmtry%clmns=clmns
    dmngmtry%rws=rws
    dmngmtry%lvls=lvls
    dmngmtry%Xd=>X                                                                                                                  !Apunta al arreglo X
    dmngmtry%Yd=>Y                                                                                                                  !Apunta al arreglo Y
    dmngmtry%Zd=>Z                                                                                                                  !Apunta al arreglo Z
    NULLIFY(dmngmtry%Xs,dmngmtry%Ys,dmngmtry%Zs)                                                                                    !Anula atributos de presicion simple
    !
END SUBROUTINE bld_domain_dp
!
!***********************************************************************************************************************************
!***********************************  DESTRUCCION DE LA VARIABLE TIPO GEOMETRIA ****************************************************
!*****************************@authores: Alvarez-Villa O.D. & Rendon-Alvarez J.P. & Perez K. ***************************************
!***********************************************************************************************************************************
!
!DESCRIPCI�N:
!Nombre: dstry_domain(destroy domain)
!Una rutina que desarma una variable tipo Domino. Cuando el Dominio se haya construido en la memoria RAM lo que hace este
!algoritmo es llenar con los valores de faltantes los atributos que corresponden a metadatos y romper los punteros de los arreglos
!de coordenadas. Es importante aclarar que los arreglos no se desalojan, y en consecuencia la liberaci�n de memoria RAM 
!no es significativa; para liberar la RAM que ocupa el Dominio deben desalojarse los objetivos a los que apuntan los arreglos de 
!coordenadas. Cuando el Dominio se haya construido en el disco duro, pueden o no romperse los v�nculos que proporcionan acceso a los 
!archivos con los arreglos de coordenadas
!
!DATOS DE ENTRADA:
!dmngmtry         : Variable tipo Dominio a destruir                                                                    {metaraster}
!switch           : Cuando es verdadera se destruye el Domino espacial, bien sea que se haya construido en la RAM           {l�gica}
!                   o en el disco duro. Cuando es falsa s�lo se destruye en la RAM
!
SUBROUTINE dstry_domain(dmngmtry,switch)    
    USE anisotype
    !
    IMPLICIT NONE
    !
    LOGICAL(LGT), INTENT(IN)    :: switch
    TYPE(domain), INTENT(INOUT) :: dmngmtry
    !
    IF (switch) THEN
        dmngmtry%exist=.false.
        dmngmtry%infile=.false.
        dmngmtry%rtdom=''
        dmngmtry%flind=-9999
        dmngmtry%rws=-9999
        dmngmtry%clmns=-9999
        dmngmtry%lvls=-9999
    END IF
    IF (ASSOCIATED(dmngmtry%Xd)) THEN
        NULLIFY(dmngmtry%Xd)
    END IF
    IF (ASSOCIATED(dmngmtry%Yd)) THEN
        NULLIFY(dmngmtry%Yd)
    END IF
    IF (ASSOCIATED(dmngmtry%Zd)) THEN
        NULLIFY(dmngmtry%Zd)
    END IF
    !
END SUBROUTINE dstry_domain
!!
!***********************************************************************************************************************************
!******************************* OBTENCION DEL CENTRO DE UN VOXEL DE UN DOMINIO ESPACIAL *******************************************
!*****************************@authores: Alvarez-Villa O.D. & Rendon-Alvarez J.P. & Perez K. ***************************************
!***********************************************************************************************************************************
!
!DESCRIPCION:
!Un conjunto de subrutinas que retornan las coordeanadas del centro de un voxel contenido en un Domino espacial, bien sea que
!este almacenado en la RAM o en un archivo (*.domnRST) en el disco duro
!
!DATOS DE ENTRADA:
!dmngmtry           : Variable tipo geometria construida en archivo                                               {geometria}
!i                  : Columna del voxel donde se encuentra el dato                                                [opcional]{entero}
!j                  : Fila del voxel donde se encuentra el dato                                                   [opcional]{entero}
!k                  : Nivel del voxel donde se encuentra el dato                                                  [opcional]{entero}
!v                  : �ndice de la variable a leer en el arreglo de datos                                         [opcional]{entero}      
!switch             : Variable logica que sirve para indicar si se desea revisar la adecuada construccion del Dominio       {l�gico}
!
!DATOS DE SALIDA:
!x,y,z              : Coordenadas leidas                                                                                      {real}
!
SUBROUTINE gt_vxlcntr_dp(dmngmtry,i,j,k,x,y,z)
    USE anisotype
    USE tpvrblchckng, ONLY : chckdmn
    USE util, ONLY : error
    !
    IMPLICIT NONE        
    !
    INTEGER(I4B), INTENT(IN)    :: i,j,k
    REAL(DP), INTENT(OUT)       :: x,y,z
    TYPE(domain), INTENT(IN)    :: dmngmtry
    !
    REAL(DP)     :: xi,yi,zi,xf,yf,zf
    !
    !Verificacion de la consistencia y pertinencia de los datos de entrada                                                                         
    IF ((i<1).OR.(i>dmngmtry%clmns)) CALL error('En gt_vxlcntr_dp: Indice del voxel esta fuera de rango')                                 
    IF ((j<1).OR.(j>dmngmtry%rws)) CALL error('En gt_vxlcntr_dp: Indice del voxel esta fuera de rango')                                         
    IF ((k<1).OR.(k>dmngmtry%lvls)) CALL error('En gt_vxlcntr_dp: Indice del voxel esta fuera de rango')                                    
    !Obtencion del valor
    xi=dmngmtry%Xd(i);xf=dmngmtry%Xd(i+1)                                                                                       
    yi=dmngmtry%Yd(j);yf=dmngmtry%Yd(j+1)                                                                                       
    zi=dmngmtry%Zd(k);zf=dmngmtry%Zd(k+1)                                                                                      
    x=0.5_dp*(xf+xi)                                                                                                          
    y=0.5_dp*(yf+yi)                                                                                                           
    z=0.5_dp*(zf+zi)           
!    WRITE(*,*)i,j,k
    !
END SUBROUTINE gt_vxlcntr_dp
!!
!***********************************************************************************************************************************
!************************************ OBTENCION DE LOS LIMITES DEL DOMINIO ESPACIAL ************************************************
!*****************************@authores: Alvarez-Villa O.D. & Rendon-Alvarez J.P. & Perez K. ***************************************
!***********************************************************************************************************************************
!
!DESCRIPCI�N:
!Un conjunto de subrutinas que retornan las coordeanadas l�mites inferiores y superiores del dominio espacial.
!
!DATOS DE ENTRADA:
!dmngmtry           : Variable tipo dominio espacial                                                                        {domain}     
!switch             : Variable l�gica que sirve para indicar si se desea revisar la adecuada construcci�n del Dominio       {l�gico}
!
!DATOS DE SALIDA:
!xi,yi,zi           : Coordenadas limItrofes inferiores del dominio espacil                                                   {real}
!xf,yf,zf           : Coordenadas limItrofes superiores del dominio espacial                                                  {real}
!
SUBROUTINE gt_dmnbnds_dp(dmngmtry,xi,yi,zi,xf,yf,zf)
    !
    USE anisotype
    USE tpvrblchckng, ONLY : chckdmn
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(OUT)       :: xi,yi,zi,xf,yf,zf
    TYPE(domain), INTENT(IN)    :: dmngmtry
    !                                                                    
    !Obtencion del valor
    xi=dmngmtry%Xd(1);yi=dmngmtry%Yd(1);zi=dmngmtry%Zd(1)                                                                      
    xf=dmngmtry%Xd(dmngmtry%clmns+1);yf=dmngmtry%Yd(dmngmtry%rws+1);zf=dmngmtry%Zd(dmngmtry%lvls+1)                             
    !IF
END SUBROUTINE gt_dmnbnds_dp
!!
!***********************************************************************************************************************************
!***************************************OBTENCION DE LA LONGITUD DE LOS LADOS DE UN VOXEL*******************************************
!*****************************@authores: Alvarez-Villa O.D. & Rendon-Alvarez J.P. & Perez K. ***************************************
!***********************************************************************************************************************************
!
!DESCRIPCION:
!Devuelve la longitud de los lados de un voxel en cada direccion cartesiana (dx,dy,dz)
!
!DATOS DE ENTRADA:
!dmngmtry           : Variable tipo geometria construida en archivo o RAM                                              {geometria}
!i                  : Columna del voxel donde se encuentra el dato                                                [opcional]{entero}
!j                  : Fila del voxel donde se encuentra el dato                                                   [opcional]{entero}
!k                  : Nivel del voxel donde se encuentra el dato                                                  [opcional]{entero}
!v                  : �ndice de la variable a leer en el arreglo de datos                                         [opcional]{entero}      
!switch             : Variable logica que sirve para indicar si se desea revisar la adecuada construccion del Dominio       {logico}
!
!DATOS DE SALIDA:
!dx,dy,dz              : Longitud de los lados del voxel i,j,k
SUBROUTINE gt_delta_dp(dmngmtry,i,j,k,dx,dy,dz)
    
    USE anisotype
    USE tpvrblchckng, ONLY : chckdmn
    USE util, ONLY : error
    !
    IMPLICIT NONE        
    !
    INTEGER(I4B), INTENT(IN)    :: i,j,k
    REAL(DP), INTENT(OUT)       :: dx,dy,dz
    TYPE(domain), INTENT(IN)   :: dmngmtry
    !
    REAL(DP)                    :: xi,yi,zi,xf,yf,zf
    !
    !Verifica que indices se encuentra dentro del dominio espacial                                                                  
    IF ((i<1).OR.(i>dmngmtry%clmns)) CALL error('En gt_vxlcntr_dp: indice de voxel fuera de rango')                                       
    IF ((j<1).OR.(j>dmngmtry%rws)) CALL error('En gt_vxlcntr_dp: indice de voxel fuera de rango')                                         
    IF ((k<1).OR.(k>dmngmtry%lvls)) CALL error('En gt_vxlcntr_dp: indice de voxel fuera de rango')                                        
    !Obtencion del valor                                                                                     
    xi=dmngmtry%Xd(i);xf=dmngmtry%Xd(i+1)                                                                                       
    yi=dmngmtry%Yd(j);yf=dmngmtry%Yd(j+1)                                                                                       
    zi=dmngmtry%Zd(k);zf=dmngmtry%Zd(k+1)                                                                                       
    dx=xf-xi                                                                                                                    
    dy=yf-yi                                                                                                                    
    dz=zf-zi                                                                                                                    
    !
END SUBROUTINE gt_delta_dp

!
!***********************************************************************************************************************************
!************************** OBTENCION DE LA DISTANCIA ENTRE LOS CENTROS DE DOS BLOQUES ADYACENTES **********************************
!*****************************@authores: Alvarez-Villa O.D. & Rendon-Alvarez J.P. & Perez K. ***************************************
!***********************************************************************************************************************************
!
!DESCRIPCION
!Permite el calculo de la distancia entre el centro de un bloque y el centro de un bloque adyacente.

!DATOS DE ENTRADA:
!dmngmtry           : Variable tipo geometria construida en archivo o RAM                                              {geometria}
!i                  : Columna del voxel donde se encuentra el dato                                                [opcional]{entero}
!j                  : Fila del voxel donde se encuentra el dato                                                   [opcional]{entero}
!k                  : Nivel del voxel donde se encuentra el dato                                                  [opcional]{entero}
!c                  : Indice que cuenta el bloque adyacenten al bloque central        
!
!DATOS DE SALIDA:
!dx,dy,dz              : Longitud de los lados del voxel i,j,k                                                                        
SUBROUTINE gt_distance_dp(dmngmtry,i,j,k,c,dist)
    USE anisotype
    USE tpvrblchckng, ONLY : chckdmn
    USE util, ONLY : error
    !
    IMPLICIT NONE        
    !
    INTEGER(I4B), INTENT(IN)    :: i,j,k
    CHARACTER(3), INTENT(IN)    :: c
    REAL(DP), INTENT(OUT)       :: dist
    TYPE(domain), INTENT(IN)    :: dmngmtry
    !
    REAL(DP)     :: xa,ya,za,xb,yb,zb,xc,yc,zc
    !
    !Verificaci�n de la consistencia y pertinencia de los datos de entrada                                                                   
    IF ((i<1).OR.(i>dmngmtry%clmns)) CALL error('in gt_vxlcntr_dp: indice del voxel esta fuera de rango')    
    IF ((j<1).OR.(j>dmngmtry%rws)) CALL error('in gt_vxlcntr_dp: indice del voxel esta fuera de rango')                                       
    IF ((k<1).OR.(k>dmngmtry%lvls)) CALL error('in gt_vxlcntr_dp: indice del voxel esta fuera de rango')                                       
    !Obtencion del valor. Verifico el bloque adyacente segun el valor que haya tomado C.
        IF (c=='i-1') THEN
            xa=dmngmtry%Xd(i-1);xb=dmngmtry%Xd(i);xc=dmngmtry%Xd(i+1)   
            dist=0.5_dp*(xb+xc)-0.5_dp*(xa+xb)
        ELSE IF (c=='i+1') THEN
            xa=dmngmtry%Xd(i);xb=dmngmtry%Xd(i+1);xc=dmngmtry%Xd(i+2)   
            dist=0.5_dp*(xb+xc)-0.5_dp*(xa+xb)            
        ELSE IF (c=='j-1') THEN
            ya=dmngmtry%Yd(j-1);yb=dmngmtry%Yd(j);yc=dmngmtry%Yd(j+1)   
            dist=0.5_dp*(yb+yc)-0.5_dp*(ya+yb)           
        ELSE IF (c=='j+1') THEN
            ya=dmngmtry%Yd(j);yb=dmngmtry%Yd(j+1);yc=dmngmtry%Yd(j+2)   
            dist=0.5_dp*(yb+yc)-0.5_dp*(ya+yb)             
        ELSE IF (c=='k-1') THEN
            za=dmngmtry%Zd(k-1);zb=dmngmtry%Zd(k);zc=dmngmtry%Zd(k+1)   
            dist=0.5_dp*(zb+zc)-0.5_dp*(za+zb)  
        ELSE IF (c=='k+1') THEN
            za=dmngmtry%Zd(k);zb=dmngmtry%Zd(k+1);zc=dmngmtry%Zd(k+2)   
            dist=0.5_dp*(zb+zc)-0.5_dp*(za+zb)   
        END IF                                                                           
    !
END SUBROUTINE gt_distance_dp
