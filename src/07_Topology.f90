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
MODULE topology 
    INTERFACE bld_topologybn
        SUBROUTINE bld_topologybn(prjctrt,dmngmtry,topo)
            USE anisotype
            IMPLICIT NONE
#include <petsc/finclude/petscsys.h>
            CHARACTER(*), INTENT(IN)          :: prjctrt   
            TYPE(domain), INTENT(IN)          :: dmngmtry
            TYPE(topology), INTENT(OUT)       :: topo
        END SUBROUTINE bld_topologybn
    END INTERFACE bld_topologybn
    
!    INTERFACE bld_topology
!        SUBROUTINE bld_topology(prjctrt,dmngmtry,topo)
!            USE anisotype
!            IMPLICIT NONE
!            CHARACTER(*), INTENT(IN)          :: prjctrt   
!            TYPE(domain), INTENT(IN)          :: dmngmtry
!            TYPE(topology), INTENT(OUT)       :: topo
!        END SUBROUTINE bld_topology
!    END INTERFACE bld_topology
END MODULE topology

!***********************************************************************************************************************************
!***********************CONSTRUCCION DE LA VARIABLE TIPO TOPOLOGIA******************************************************************
!***************************@author Álvarez-Villa O.D. & Perez K.*******************************************************************
!***********************************************************************************************************************************    
!
!DESCRIPCION:
!NOMBRE: bld_topology 
!Subrutina que llena los atibutos de la variable tipo toplogia a partir de informacion almacenada en numeros binarios.
!
!DATOS DE ENTRADA:
!prjctrt        :Ruta del proyecto de trabajo
!dmngmtry       :Variable tipo geometria
!
!DATOS DE SALIDA:
!topo           :Variable tipo de topologia cosntruida.

!OBSERVACIONES: Se crean primero los arreglos en RAM que contienen el numero de celdas de cada tipo y los arreglos con las 
!coordenadas

SUBROUTINE bld_topologybn(prjctrt,dmngmtry,topo)
    USE anisotype
    USE iodata, ONLY : gtvlr_bn
    !
    IMPLICIT NONE
#include <petsc/finclude/petscsys.h>
    PetscErrorCode                              :: ierr
    PetscMPIInt                                 :: tprocess, iprocess,status(MPI_STATUS_SIZE)
    !
    CHARACTER(*), INTENT(IN)                    :: prjctrt   
    TYPE(domain), INTENT(IN)                    :: dmngmtry
    TYPE(topology), INTENT(INOUT)               :: topo
    !   
    CHARACTER(200)                              :: rtact,rtcth,rtdph,rtk,rtch,rteh
    INTEGER(I4B)                                :: actk,dphk,cont,cero,act,cth,dph 
    INTEGER(I4B)                                :: uact,udph,utop,uch,ueh,i,j,k,a
    !
    PARAMETER(uact=30); PARAMETER(udph=50)
    PARAMETER(utop=55); PARAMETER(uch=57); PARAMETER(ueh=58)
    !
    cero=0; cont=0;
    !
    CALL MPI_comm_rank(MPI_COMM_WORLD,iprocess,ierr)
    CALL MPI_comm_size(MPI_COMM_WORLD,tprocess,ierr)
    ! Barrier: Sincroniza todos los procesadores al inicio, esto para ubicar más fácil posibles errores.
    CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
    IF (iprocess==0) THEN
        WRITE (*,*) '*---|Inicia construccion de la topologia|---*'
        !
        rtact = trim(prjctrt)//'/sanpck.act'   !Archivo que contiene los bloques activos del modelo
        rtdph = trim(prjctrt)//'/sanpck.dph'   !Archivo que contiene los bloques con CC de nivel externo
        rtk   = trim(prjctrt)//'/sanpck.top'   !Archivo que contiene indices topologicos para bloques activos
        rtch  = trim(prjctrt)//'/sanpck.toc'   !Archivo que contiene indices topologicos para CC de nivel impuesto
        rteh  = trim(prjctrt)//'/sanpck.toe'   !Archivo que contiene indices topologicos para CC de nivel externo
        OPEN(uact, FILE=rtact, STATUS='OLD', ACCESS='DIRECT', FORM='UNFORMATTED', RECL=4*6)
        OPEN(udph, FILE=rtdph, STATUS='OLD', ACCESS='DIRECT', FORM='UNFORMATTED', RECL=4*6)
        OPEN(utop, FILE=rtk, STATUS='REPLACE', ACCESS='DIRECT', FORM='UNFORMATTED', RECL=10*4)
        OPEN(uch,  FILE=rtch, STATUS='REPLACE', ACCESS='DIRECT', FORM='UNFORMATTED', RECL=10*4)
        OPEN(ueh,  FILE=rteh, STATUS='REPLACE', ACCESS='DIRECT', FORM='UNFORMATTED', RECL=10*4)
        !
        act=0; cth=0; dph=0 
        !NIveles paralelos a Z
        DO k=1,dmngmtry%lvls
            !Filas paralelas a Y
            DO j=1,dmngmtry%rws
                !Columnas paralelas a X
                DO i=1,dmngmtry%clmns
                    !
                    cont=cont+1
                    CALL gtvlr_bn(uact,dmngmtry,i,j,k,actk)
                    CALL gtvlr_bn(udph,dmngmtry,i,j,k,dphk)    
                    !
                    !Impresion de los archivos de indices topologicos tanto para
                    !los bloques activos como para las condiciones de frontera
                    IF (actk==1) THEN
                        act=act+1                        
                        WRITE(utop,REC=cont)act
                        WRITE(uch,REC=cont)cero
                        IF(dphk==3) THEN
                            dph=dph+1
                            WRITE(ueh,REC=cont)dph
                        ELSE
                            WRITE(ueh,REC=cont)cero
                        END IF 
                    ELSE IF (actk==2) THEN
                        cth=cth+1
                        WRITE(utop,REC=cont)cero
                        WRITE(ueh,REC=cont)cero
                        WRITE(uch,REC=cont)cth
                    ELSE
                        WRITE(utop,REC=cont)cero
                        WRITE(uch,REC=cont)cero
                        WRITE(ueh,REC=cont)cero
                    END IF 
                    !
                END DO 
            END DO 
        END DO 
        CLOSE(utop);CLOSE(uch); CLOSE(ueh)
    END IF

    !
    !Le doy tamaño a los arreglos topol�gicos
    ALLOCATE(topo%ijactv(act,3))    !El arreglo de bloques activos
    ALLOCATE(topo%ijcth(cth,3))    !El arreglo de bloques de nivel prescrito
    ALLOCATE(topo%ijdph(dph,3))    !El arreglo de bloques de nivel externo
    act=0; cth=0; dph=0      !Inicializamiento nulo de vectores
    !
    ! Todos los procesadores tienen la información de la topología globalmente
    !Se recorre la geometr�a del acuifero y se llenan los arreglos topol�gicos
    DO k=1,dmngmtry%lvls
        DO j=1,dmngmtry%rws
            DO i=1,dmngmtry%clmns
                IF (iprocess==0) THEN
                    CALL gtvlr_bn(uact,dmngmtry,i,j,k,actk)
                    CALL gtvlr_bn(udph,dmngmtry,i,j,k,dphk)
                    DO a=1,tprocess-1
                        CALL MPI_SEND(actk,1,MPI_INTEGER,a,300,MPI_COMM_WORLD,ierr)
                        CALL MPI_SEND(dphk,1,MPI_INTEGER,a,302,MPI_COMM_WORLD,ierr)
                    END DO
                ELSE
                    CALL MPI_RECV(actk,1,MPI_INTEGER,MPI_ANY_SOURCE,300,MPI_COMM_WORLD, status, ierr)
                    CALL MPI_RECV(dphk,1,MPI_INTEGER,MPI_ANY_SOURCE,302,MPI_COMM_WORLD, status, ierr)
                END IF
                !
                IF (actk==1)THEN
                    act=act+1
                    topo%ijactv(act,1)=i
                    topo%ijactv(act,2)=j
                    topo%ijactv(act,3)=k
                ELSE IF (actk==2) THEN
                    cth=cth+1
                    topo%ijcth(cth,1)=i
                    topo%ijcth(cth,2)=j
                    topo%ijcth(cth,3)=k
                END IF 
                !
                IF(dphk==3) THEN
                    dph=dph+1
                    topo%ijdph(dph,1)=i
                    topo%ijdph(dph,2)=j
                    topo%ijdph(dph,3)=k
                END IF 
                !
            END DO 
        END DO 
    END DO 
    !Se asignan los apuntadores de la variable tipo topology a los target credos en las subrutinas anteriores.
    topo%act=act
    topo%cth=cth
    topo%dph=dph
    topo%hastopo=.true.
    !
    !
    ! Barrier: Sincroniza todos los procesadores al final, esto para ubicar más fácil posibles errores.
    CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
    IF (iprocess==0) THEN
        CLOSE(uact); CLOSE(udph);
        WRITE (*,*) '*---|Finaliza construccion de la topologia|---*'
    END IF
    !
END SUBROUTINE  bld_topologybn

!DESCRIPCION:
!NOMBRE: bld_topologybn 
!Subrutina que llena los atibutos de la variable tipo toplogia a partir de informacion alamacenada en archivos de texto.
!
!DATOS DE ENTRADA:
!prjctrt        :Ruta del proyecto de trabajo
!cfnd           :La geometria de un acuifero confinado lineal.
!
!DATOS DE SALIDA:
!topo           :Variable tipo de topologia cosntruida.
!
!SUBROUTINE bld_topology(prjctrt,dmngmtry,topo)
!    USE anisotype
!    USE iodata, ONLY : gtvlr_bn
!    !
!    IMPLICIT NONE
!    !
!    CHARACTER(*), INTENT(IN)                                     :: prjctrt   
!    TYPE(domain), INTENT(IN)                                     :: dmngmtry
!    TYPE(topology), INTENT(OUT)                                  :: topo
!    !   
!    CHARACTER(200)                                               :: rtact,rtcth,rtdph,rtk,rtch,rteh
!    INTEGER(I4B)                                                 :: actk,dphk,cont,cero,act,cth,dph 
!    INTEGER(I4B)                                                 :: uact,udph,utop,uch,ueh,i,j,k
!    INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE, TARGET            :: actv,cthv,dphv
!    !
!    PARAMETER(uact=30); PARAMETER(udph=50)
!    PARAMETER(utop=55); PARAMETER(uch=57); PARAMETER(ueh=58)
!    !
!    cero=0; cont=0;
!    !
!    WRITE (*,*) '*---|Inicia construccion de la topologia|---*'
!    !
!    rtact = trim(prjctrt)//'/sanpck.act'   !Archivo que contiene los bloques activos del modelo
!    rtdph = trim(prjctrt)//'/sanpck.dph'   !Archivo que contiene los bloques con CC de nivel externo
!    rtk   = trim(prjctrt)//'/sanpck.top'   !Archivo que contiene indices topologicos para bloques activos
!    rtch  = trim(prjctrt)//'/sanpck.toc'   !Archivo que contiene indices topologicos para CC de nivel impuesto
!    rteh  = trim(prjctrt)//'/sanpck.toe'   !Archivo que contiene indices topologicos para CC de nivel externo
!    OPEN(uact, FILE=rtact, STATUS='OLD', ACCESS='DIRECT', RECL=6, FORM='FORMATTED')
!    OPEN(udph, FILE=rtdph, STATUS='OLD', ACCESS='DIRECT', RECL=6, FORM='FORMATTED')
!    OPEN(utop, FILE=rtk, STATUS='REPLACE', ACCESS='DIRECT', RECL=11, FORM='FORMATTED')
!    OPEN(uch,  FILE=rtch, STATUS='REPLACE', ACCESS='DIRECT', RECL=11, FORM='FORMATTED')
!    OPEN(ueh,  FILE=rteh, STATUS='REPLACE', ACCESS='DIRECT', RECL=11, FORM='FORMATTED')
!    !
!    act=0; cth=0; dph=0 
!    !NIveles paralelos a Z
!    DO k=1,dmngmtry%lvls
!        !Filas paralelas a Y
!        DO j=1,dmngmtry%rws
!            !Columnas paralelas a X
!            DO i=1,dmngmtry%clmns
!                !
!                cont=cont+1
!                CALL gtvlr(uact,dmngmtry,i,j,k,actk)
!                CALL gtvlr(udph,dmngmtry,i,j,k,dphk)    
!                !
!                !Impresion de los archivos de indices topologicos tanto para
!                !los bloques activos como para las CC de nivel impuesto
!                IF (actk==1) THEN
!                    act=act+1                        
!                    WRITE(utop,'((I10),(A1))',REC=cont)act,13
!                    WRITE(uch,'((I10),(A1))',REC=cont)cero,13
!                    IF(dphk==3) THEN
!                        dph=dph+1
!                        WRITE(ueh,'((I10),(A1))',REC=cont)dph,13
!                    ELSE
!                        WRITE(ueh,'((I10),(A1))',REC=cont)cero,13
!                    END IF 
!                ELSE IF (actk==2) THEN
!                    cth=cth+1
!                    WRITE(utop,'((I10),(A1))',REC=cont)cero,13
!                    WRITE(ueh,'((I10),(A1))',REC=cont)cero,13
!                    WRITE(uch,'((I10),(A1))',REC=cont)cth,13
!                ELSE
!                    WRITE(utop,'((I10),(A1))',REC=cont)cero,13
!                    WRITE(uch,'((I10),(A1))',REC=cont)cero,13
!                    WRITE(ueh,'((I10),(A1))',REC=cont)cero,13
!                END IF 
!                !
!            END DO 
!        END DO 
!    END DO 
!    CLOSE(utop);CLOSE(uch); CLOSE(ueh)
!    !
!    !Le doy tama�o a los arreglos topol�gicos
!    ALLOCATE(actv(act,3))   !El arreglo de bloques activos
!    ALLOCATE(cthv(cth,3))    !El arreglo de bloques de nivel prescrito
!    ALLOCATE(dphv(dph,3))    !El arreglo de bloques de nivel externo
!    act=0; cth=0; dph=0  !Inicializamiento nulo de vectores
!    !
!    !Se recorre la geometr�a del acuifero y se llenan los arreglos topol�gicos
!    DO k=1,dmngmtry%lvls
!        DO j=1,dmngmtry%rws
!            DO i=1,dmngmtry%clmns
!                !
!                CALL gtvlr(uact,dmngmtry,i,j,k,actk)
!                CALL gtvlr(udph,dmngmtry,i,j,k,dphk)
!                !
!                IF (actk==1)THEN
!                    act=act+1
!                    actv(act,1)=i
!                    actv(act,2)=j
!                    actv(act,3)=k
!                ELSE IF (actk==2) THEN
!                    cth=cth+1
!                    cthv(cth,1)=i
!                    cthv(cth,2)=j
!                    cthv(cth,3)=k
!                END IF 
!                !
!                IF(dphk==3) THEN
!                    dph=dph+1
!                    dphv(dph,1)=i
!                    dphv(dph,2)=j
!                    dphv(dph,3)=k
!                END IF 
!                !
!            END DO 
!        END DO 
!    END DO 
!    !Se asignan los apuntadores de la variable tipo topology a los target credos en las subrutinas anteriores.
!    topo%act=act
!    topo%cth=cth
!    topo%dph=dph
!    topo%hastopo=.true.
!    topo%ijactv=>actv
!    topo%ijcth=>cthv
!    topo%ijdph=>dphv
!    !
!    CLOSE(uact); CLOSE(udph);
!    WRITE (*,*) '*---|Finaliza construccion de la topologia|---*'
!    !
!END SUBROUTINE  bld_topology
!
