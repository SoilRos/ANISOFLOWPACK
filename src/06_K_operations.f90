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

MODULE k_operations
    INTERFACE K_tensor
        SUBROUTINE K_tensor_bn(prjctrt,dmngmtry,kxx,kyy,kzz,kxy,kxz,kyz)
            USE anisotype
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN)          :: prjctrt   
            TYPE(domain), INTENT(IN)          :: dmngmtry
            REAL(DP), DIMENSION(:),ALLOCATABLE, INTENT(OUT)                :: kxx,kyy,kzz,kxy,kxz,kyz
        END SUBROUTINE K_tensor_bn
    END INTERFACE K_tensor
END MODULE k_operations

SUBROUTINE K_tensor_bn(prjctrt,dmngmtry,kxx,kyy,kzz,kxy,kxz,kyz)
    USE anisotype
    USE geometry
    USE iodata,     ONLY : gtvlr_bn
    USE util
    !
    IMPLICIT NONE
    !
    CHARACTER(*), INTENT(IN)                                       :: prjctrt
    TYPE(domain), INTENT(IN)                                       :: dmngmtry
    REAL(DP), DIMENSION(:),ALLOCATABLE, INTENT(OUT)                :: kxx,kyy,kzz,kxy,kxz,kyz
    !
    CHARACTER(200)                                                 :: rttrns,rtalf,rtbet,rttet,rtkx,rtky,rtkz,rtxy,rtxz,rtyz
    REAL(DP)                                                       :: ka,kb,kc,kd,ke,kf,kg,kh,ki,al,be,te
    INTEGER(I4B)                                                   :: u,ua,ub,ut,uc,ud,ue,uf,ug,uh,i,j,k,n
    INTEGER(I4B)                                                   :: alfa,beta,teta,cont
    !
    PARAMETER(u=37,ua=65,ub=58,ut=33,uc=38,ud=39,ue=45,uf=74,ug=21,uh=22)
    !
    WRITE (*,*) '*---|Inicia calculo de tensor anisotropo de conductividades |---*'
    !Propiedades para determinar las direcciones de anisotropia
    rttrns = trim(prjctrt)//'/sanpck.cvt'  !Archivo de conductividad de bloque
    rtalf  = trim(prjctrt)//'/sanpck.alf'  !Archivo de longitudes del rio
    rtbet  = trim(prjctrt)//'/sanpck.bet'  !Archivo de anchos del rio
    rttet  = trim(prjctrt)//'/sanpck.tet'  !Archivo de espesor en el lecho de rio
    OPEN(u, FILE=rttrns, STATUS='OLD', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
    OPEN(ua, FILE=rtalf, STATUS='OLD', ACCESS='DIRECT', RECL=4*6, FORM='UNFORMATTED')
    OPEN(ub, FILE=rtbet, STATUS='OLD', ACCESS='DIRECT', RECL=4*6, FORM='UNFORMATTED')
    OPEN(ut, FILE=rttet, STATUS='OLD', ACCESS='DIRECT', RECL=4*6, FORM='UNFORMATTED')
    !
!    rtkx = trim(prjctrt)//'/sanpck.kxx' !Archivo con conductividades en direccion xx
!    rtky = trim(prjctrt)//'/sanpck.kyy' !Archivo con conductividades en direccion yy
!    rtkz = trim(prjctrt)//'/sanpck.kzz' !Archivo con conductividades en direccion zz
!    rtxy= trim(prjctrt)//'/sanpck.kxy' !Archivo con conductividades en direccion xy
!    rtxz = trim(prjctrt)//'/sanpck.kxz' !Archivo con conductividades en direccion xz
!    rtyz = trim(prjctrt)//'/sanpck.kyz' !Archivo con conductividades en direccion yz
!    OPEN(uc, FILE=rtkx, STATUS='REPLACE', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
!    OPEN(ud, FILE=rtky, STATUS='REPLACE', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
!    OPEN(ue, FILE=rtkz, STATUS='REPLACE', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
!    OPEN(uf, FILE=rtxy, STATUS='REPLACE', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
!    OPEN(ug, FILE=rtxz, STATUS='REPLACE', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
!    OPEN(uh, FILE=rtyz, STATUS='REPLACE', ACCESS='DIRECT', RECL=8*18, FORM='UNFORMATTED')
    cont=0
    n=dmngmtry%clmns*dmngmtry%rws*dmngmtry%lvls
    ALLOCATE(kxx(n));                         kxx=0.0_dp
    ALLOCATE(kyy(n));                         kyy=0.0_dp
    ALLOCATE(kzz(n));                         kzz=0.0_dp
    ALLOCATE(kxy(n));                         kxy=0.0_dp
    ALLOCATE(kxz(n));                         kxz=0.0_dp
    ALLOCATE(kyz(n));                         kyz=0.0_dp
    
    DO k=1,dmngmtry%lvls
        DO j=1,dmngmtry%rws
            DO i=1,dmngmtry%clmns
                cont=cont+1
                Ka=0.0_dp; Kb=0.0_dp; Kc=0.0_dp; Kd=0.0_dp; Ke=0.0_dp; Kf=0.0_dp; Kg=0.0_dp
!                kxx=0.0_dp; kyy=0.0_dp; kzz=0.0_dp; kxy=0.0_dp; kxz=0.0_dp; kyz=0.0_dp
                alfa=0; beta=0; teta=0
                !Calculo conductividades de bloques adyacentes al bloque i,j,k
                CALL gtvlr_bn(u,dmngmtry,i,j,k,ka) 
                !Calculo de conductividades sobre ejes cartesianos en direccion X
                IF ( i==dmngmtry%clmns ) THEN
                    CALL gtvlr_bn(u,dmngmtry,i-1,j,k,kb)
                    kxx(cont)=armncmn(ka,kb)
                ELSE  
                    CALL gtvlr_bn(u,dmngmtry,i+1,j,k,kc)
                    kxx(cont)=armncmn(ka,kc)
                END IF
                !Calculo de conductividades sobre ejes cartesianos en direccion Y
                IF ( j==dmngmtry%rws ) THEN
                    CALL gtvlr_bn(u,dmngmtry,i,j-1,k,kd)
                    kyy(cont)=armncmn(ka,kd)
                ELSE 
                    CALL gtvlr_bn(u,dmngmtry,i,j+1,k,ke)
                    kyy(cont)=armncmn(ka,ke)
                END IF
                !Calculo de conductividades sobre ejes cartesianos en direccion Z
                IF ( k==dmngmtry%lvls )  THEN
                    CALL gtvlr_bn(u,dmngmtry,i,j,k-1,kf)
                    kzz(cont)=armncmn(ka,kf)
                ELSE 
                    CALL gtvlr_bn(u,dmngmtry,i,j,k+1,kg)
                    kzz(cont)=armncmn(ka,kg)
                END IF
!                WRITE(*,*)kxx,kyy,kzz,'kxx,kyy,kzz'
                !Calculo de conductividades de los ejes principales y de las direcciones paralelas a planos (cortantes)
                CALL gtvlr_bn(ua,dmngmtry,i,j,k,alfa) !Rotacion de eje X
                CALL gtvlr_bn(ub,dmngmtry,i,j,k,beta) !Rotacion de eje Y
                CALL gtvlr_bn(ut,dmngmtry,i,j,k,teta) !Rotacion de eje Z
!                WRITE(*,*)alfa,beta,teta,'alfa,teta,beta'
                !Convierto angulos a radianes
                al=PI_D*real(alfa)/180 
                be=PI_D*real(beta)/180
                te=PI_D*real(teta)/180
                kh=-kxx(cont)/(-cos(te)*cos(be)+sin(te)*cos(be)-sin(be)) !Direccion principal de la conductividad X
                ki=(-kyy(cont)-kh*sin(te)*cos(be))/(cos(be)*sin(al)-sin(te)*cos(al)+cos(te)*sin(be)*sin(al))!Direccion principal de la conductividad Y
                kxy(cont)=abs((-sin(te)*cos(al)+cos(te)*sin(be)*sin(al))*ki)
                kxz(cont)=abs(-sin(be)*kh)
                kyz(cont)=abs(cos(be)*sin(al)*ki)
!                WRITE(*,*)kxy,kxz,kyz,'kxy,kxz,kyz'
!                WRITE(*,*)kh,ki
!                WRITE(uc, rec=cont)kxx  !Escribe archivo de conductividad xx
!                WRITE(ud, rec=cont)kyy  !Escribe archivo de conductividad yy
!                WRITE(ue, rec=cont)kzz  !Escribe archivo de conductividad zz
!                WRITE(uf, rec=cont)kxy  !Escribe archivo de conductividad xy
!                WRITE(ug, rec=cont)kxz  !Escribe archivo de conductividad xz
!                WRITE(uh, rec=cont)kyz  !Escribe archivo de conductividad yz
                !
            END DO
        END DO
    END DO
    !
    CLOSE(u)
    CLOSE(ua) 
    CLOSE(ub) 
    CLOSE(ut) 
!    CLOSE(uc) 
!    CLOSE(ud) 
!    CLOSE(ue) 
!    CLOSE(uf) 
!    CLOSE(ug) 
!    CLOSE(uh) 
    !
    WRITE (*,*) '*---|Finaliza calculo de tensor anisotropo de conductividades |---*'
    !
END SUBROUTINE K_tensor_bn