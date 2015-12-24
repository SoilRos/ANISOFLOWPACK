!!-----------------------------------------------------------------------------------------------------------------------------------
!!Copyright (c) 2013, Alvarez-Villa O.D., Perez Tevin  GOTTA Ingenieria SAS, All Right Reserved.
!!
!!This Program is distributed in the hope that they will be used but WITHOUT ANY WARRANTY. No autor or distributor accepts 
!!responsability to anyone for consequences of using them or for whether they serve !any particular porpouse or work at all, unless he
!!says so in writing. Everyone is granted permission to copy, modify an distribute this program, but only under the conditions that 
!!this notice and above copyright notice remains intact.
!!-----------------------------------------------------------------------------------------------------------------------------------
!!-----------------------------------------------------------------------------------------------------------------------------------
!!Este modulo contiene las variables tipo usadas en ANISOFLOWPACK. Las variables tipo implementadas han sido adaptadas partir de los
!!modulos de variables tipo de los siguientes trabajos:
!!1. Festellustype del software de interpolacion espacial usando tecnicas geoestadisticas Festellus. 2014, Rendon-Alvarez J.P.
!!2. Eigentype del software para modelar flujo subterraneo FDPACK. 2013, Alvarez-Villa O.D.
!!-----------------------------------------------------------------------------------------------------------------------------------
!
!MODULE linsolver
!        !!
!    INTERFACE lubksb
!        SUBROUTINE lubksb_sp(a,indx,b)
!            USE anisotype
!            IMPLICIT NONE
!            REAL(SP), DIMENSION(:,:), INTENT(IN)   :: a
!            INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
!            REAL(SP), DIMENSION(:), INTENT(INOUT)  :: b
!        END SUBROUTINE lubksb_sp
!        !BL
!        SUBROUTINE lubksb_dp(a,indx,b)
!            USE anisotype
!            IMPLICIT NONE
!            REAL(DP), DIMENSION(:,:), INTENT(IN)   :: a
!            INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
!            REAL(DP), DIMENSION(:), INTENT(INOUT)   :: b
!        END SUBROUTINE lubksb_dp
!    END INTERFACE
!	!!
!    INTERFACE ludcmp
!        SUBROUTINE ludcmp_sp(a,indx,d)
!            USE anisotype
!            IMPLICIT NONE
!            REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
!            INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
!            REAL(SP), INTENT(OUT)                   :: d
!        END SUBROUTINE ludcmp_sp
!        !BL
!        SUBROUTINE ludcmp_dp(a,indx,d)
!            USE anisotype
!            IMPLICIT NONE
!            REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
!            INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
!            REAL(DP), INTENT(OUT)                   :: d
!        END SUBROUTINE ludcmp_dp		
!    END INTERFACE
!
!END MODULE linsolver
!
!!
!!***********************************************************************************************************************************
!!******************DESCOMPOSICION LU y SUSTITUCION HACIA ATRAS PARA AX=b************************************************************
!!****************@authores: Press et al. (1986) @modified: Alvarez-Villa O.D.*******************************************************
!!***********************************************************************************************************************************
!!
!!Dada una matriz de entrada de tamaño NxN, esta rutina la reemplaza por su descomposicion LU de una permutacion por filas de 
!!ella misma. Los parametros son:
!!a          :Matriz a la cual se  le realiza la descomposicion
!!Salidas:
!!a          :Descomposicion LU de a (entrada)
!!indx       :Es un vector de longitud N que contiene las permutaciones por fila realizadas
!!d          :Es + o - 1 dependiendo de si el intercambio es par o impar
!!
!SUBROUTINE ludcmp_sp(a,indx,d)
!    USE anisotype
!    USE util, ONLY : assert_eq,imaxloc,error,outerprod,swap
!    !
!    IMPLICIT NONE
!    !
!    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
!    INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
!    REAL(SP), INTENT(OUT)                   :: d   
!    !
!    REAL(SP), DIMENSION(SIZE(a,1)) :: vv                                               !vv almacena el escalamiento implicito... 
!    REAL(SP), PARAMETER            :: TINY=1.0e-20_sp                                  !...de cada columna
!    INTEGER(I4B)                   :: j,n,imax
!    !
!    n=assert_eq(SIZE(a,1),SIZE(a,2),SIZE(indx),'ludcmp')
!    d=1.0                                                                              !No hay intercambios de fila hasta ahora
!    vv=MAXVAL(ABS(a),dim=2)                                                            !Loop sobre las filas para obtener la... 
!    IF (ANY(vv == 0.0)) THEN                                                           !...informacion de escalamiento implicito                      
!        CALL error('in ludcmp: singular matrix')                                       !Si hay una fila con ceros, error
!    END IF    
!    vv=1.0_sp/vv                                                                       !Calcula el escalamiento
!    !
!    LOOPPPAL: DO j=1,n                                                                 !Loop principal
!        !
!        imax=(j-1)+imaxloc(vv(j:n)*ABS(a(j:n,j)))                                      !Buscar la fila pivote
!        IF (j /= imax) THEN                                                            !Intercambio de filas, de ser necesario
!            CALL swap(a(imax,:),a(j,:))
!            d=-d
!            vv(imax)=vv(j)
!        END IF
!        indx(j)=imax
!        IF (a(j,j) == 0.0) a(j,j)=TINY                                                 !Si el pivote es cero es una matriz singular	    
!        a(j+1:n,j)=a(j+1:n,j)/a(j,j)                                                   !Division por el pivote	    
!        a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))                 !Reducir la submatriz que queda
!        !                                         
!    END DO LOOPPPAL
!    !
!END SUBROUTINE ludcmp_sp
!!BL
!SUBROUTINE ludcmp_dp(a,indx,d)
!    USE anisotype
!    USE util, ONLY : assert_eq,imaxloc,error,outerprod,swap
!    !
!    IMPLICIT NONE
!    !
!    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
!    INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
!    REAL(DP), INTENT(OUT)                   :: d    
!    !
!    REAL(DP), DIMENSION(SIZE(a,1)) :: vv                                                !vv almacena el escalamiento implicito...     
!    REAL(DP), PARAMETER            :: TINY=1.0e-20_dp                                   !...de cada columna
!    INTEGER(I4B)                   :: j,n,imax
!    !
!    n=assert_eq(SIZE(a,1),SIZE(a,2),SIZE(indx),'ludcmp')
!    d=1.0                                                                               !No hay intercambios de fila hasta ahora
!    vv=MAXVAL(ABS(a),dim=2)                                                             !Loop sobre las filas para obtener la.. 
!    IF (ANY(vv == 0.0)) THEN                                                            !...informacion de escalamiento implicito                      
!        CALL error('in ludcmp: singular matrix')                                        !Si hay una fila con ceros, error
!    END IF    
!    vv=1.0_sp/vv                                                                        !Calcula el escalamiento
!    !
!    LOOPPPAL: DO j=1,n                                                                  !Loop principal
!        !
!        imax=(j-1)+imaxloc(vv(j:n)*ABS(a(j:n,j)))                                       !Buscar la fila pivote
!        IF (j /= imax) THEN                                                             !Intercambio de filas, de ser necesario
!            CALL swap(a(imax,:),a(j,:))
!            d=-d
!            vv(imax)=vv(j)
!        END IF
!        indx(j)=imax
!        IF (a(j,j) == 0.0) a(j,j)=TINY                                                  !Si el pivote es cero es una matriz singular	    
!        a(j+1:n,j)=a(j+1:n,j)/a(j,j)                                                    !Division por el pivote	    
!        a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))                  !Reducir la submatriz que queda
!        !                                         
!    END DO LOOPPPAL
!    !
!END SUBROUTINE ludcmp_dp
!!
!!Resuelve un sistema lineal de la forma AX=b mediante sustitucion hacia atras. Los parametros son:
!!a          :Una matriz en formato descompisicion LU, para proceder a sustitucion
!!indx       :Vector de permutaciones estimado mediante ludcmp
!!b          :Vector de terminos independientes, de longitud N
!!Salida
!!b          :Los componentes de b en la salida son la soluci�n X del sistema.
!!
!SUBROUTINE lubksb_sp(a,indx,b)
!    USE anisotype
!    USE util, ONLY : assert_eq
!    !
!    IMPLICIT NONE
!    REAL(SP), DIMENSION(:,:), INTENT(IN)   :: a
!    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
!    REAL(SP), DIMENSION(:), INTENT(INOUT)  :: b
!    !
!    INTEGER(I4B) :: i,n,ii,ll
!    REAL(SP)     :: summ
!    !
!    n=assert_eq(SIZE(a,1),SIZE(a,2),SIZE(indx),'lubksb')
!    ii=0                                                                                !Cuando ii es positivo se convierte...
!    !                                                                                   !...en el indice de movimiento de b
!    LOOPPPAL: DO i=1,n
!        ll=indx(i)
!        summ=b(ll)
!        b(ll)=b(i)
!        IF (ii /= 0) THEN
!            summ=summ-DOT_PRODUCT(a(i,ii:i-1),b(ii:i-1))
!        ELSE IF (summ /= 0.0) THEN	        
!            ii=i                                                                       !Se encontro un elemento no nulo,... 
!        END IF                                                                         !...entonces se realiza producto punto
!        b(i)=summ
!    END DO LOOPPPAL
!    DO i=n,1,-1                                                                        !La sustitucion hacia atras        
!        b(i) = (b(i)-DOT_PRODUCT(a(i,i+1:n),b(i+1:n)))/a(i,i)
!    END DO
!    !
!END SUBROUTINE lubksb_sp
!!BL
!SUBROUTINE lubksb_dp(a,indx,b)
!    USE anisotype
!    USE util, ONLY : assert_eq
!    !
!    IMPLICIT NONE
!    REAL(DP), DIMENSION(:,:), INTENT(IN)   :: a
!    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
!    REAL(DP), DIMENSION(:), INTENT(INOUT)  :: b
!    !
!    INTEGER(I4B) :: i,n,ii,ll
!    REAL(DP) :: summ
!    !
!    n=assert_eq(SIZE(a,1),SIZE(a,2),SIZE(indx),'lubksb')
!    ii=0                                                                               !Cuando ii es positivo se convierte
!    !                                                                                  !en el indice de movimiento de b
!    LOOPPPAL: DO i=1,n
!        ll=indx(i)
!        summ=b(ll)
!        b(ll)=b(i)
!        IF (ii /= 0) THEN
!            summ=summ-DOT_PRODUCT(a(i,ii:i-1),b(ii:i-1))
!        ELSE IF (summ /= 0.0) THEN	        
!            ii=i                                                                       !Se encontro un elemento no nulo, 
!        END IF                                                                         !entonces se realiza producto punto
!        b(i)=summ
!    END DO LOOPPPAL
!    DO i=n,1,-1                                                                        !La sustitucion hacia atras        
!        b(i) = (b(i)-DOT_PRODUCT(a(i,i+1:n),b(i+1:n)))/a(i,i)
!    END DO
!    !
!END SUBROUTINE lubksb_dp
!!