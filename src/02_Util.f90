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
MODULE util
    USE anisotype
    !
    IMPLICIT NONE
    
    !*******************************************************************************************************************************
    !*******************************************************************************************************************************
    !RUTINAS QUE REVISAN ARGUMENTOS Y MANEJAN LOS ERRORES
    !*******************************************************************************************************************************
    !*******************************************************************************************************************************
    !
    !El programa embebido muere con un mensaje de error si uno de los argumentos logicos es falso.
    !
    INTERFACE assert
            MODULE PROCEDURE assert1,assert2,assert3,assert4,assert_v
    END INTERFACE
    !
    !El programa embebido muere con un mensaje de error si uno de los argumentos enteros no son iguales al primero. De otra forma, 
    !devuelve el primer argumento.
    !
    INTERFACE assert_eq
            MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
    END INTERFACE	
    !
    !Rellena el almacenamiento de un nuevo arreglo del tamaño especificado por
    !los enteros n, m, .....(igual al tamaño del puntero p). Luego copia el
    !contenido del puntero p en el nuevo almacenamiento.
    !
    INTERFACE reallocate
		MODULE PROCEDURE reallocate_rv,reallocate_dv,reallocate_rm,&
			reallocate_dm,reallocate_iv,reallocate_im,reallocate_hv
    END INTERFACE
    !
    CONTAINS
    !
    !*******************************************************************************************************************************
    !RUTINAS QUE REVISAN ARGUMENTOS Y ARROJAN ERRORES
    !*******************************************************************************************************************************
    !
    !ASSERTS
    !
    SUBROUTINE assert1(n1,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN)          :: n1
        IF (.not. n1) THEN
            WRITE (*,*) 'error: an assertion failed with this tag:',string
            STOP 'program terminated by assert1'
        END IF
    END SUBROUTINE assert1

    SUBROUTINE assert2(n1,n2,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN)          :: n1,n2
        IF (.not. (n1 .and. n2)) THEN
            WRITE (*,*) 'error: an assertion failed with this tag:',string
            STOP 'program terminated by assert2'
        END IF
    END SUBROUTINE assert2

    SUBROUTINE assert3(n1,n2,n3,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN)          :: n1,n2,n3
        IF (.not. (n1 .and. n2 .and. n3)) THEN
            WRITE (*,*) 'error: an assertion failed with this tag:',string
            STOP 'program terminated by assert3'
        END IF
    END SUBROUTINE assert3

    SUBROUTINE assert4(n1,n2,n3,n4,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN)          :: n1,n2,n3,n4
        IF (.not. (n1 .and. n2 .and. n3 .and. n4)) THEN
            WRITE (*,*) 'error: an assertion failed with this tag:',string
            STOP 'program terminated by assert4'
        END IF
    END SUBROUTINE assert4

    SUBROUTINE assert_v(n,string)
        CHARACTER(LEN=*), INTENT(IN)      :: string
        LOGICAL, DIMENSION(:), INTENT(IN) :: n
        IF (.not. all(n)) THEN
            WRITE (*,*) 'error: an assertion failed with this tag:',string
            STOP 'program terminated by assert_v'
        END IF
    END SUBROUTINE assert_v
    !
    !ASSERT_EQs
    !
    FUNCTION assert_eq2(n1,n2,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        INTEGER, INTENT(IN)          :: n1,n2
        INTEGER                      :: assert_eq2
        IF (n1 == n2) THEN
            assert_eq2=n1
        ELSE
            WRITE (*,*) 'error: an assert_eq failed with this tag:',string
            WRITE (*,*) 'push any key to finish'
            READ  (*,*)
            STOP 'program terminated by assert_eq2'
        END IF
    END FUNCTION assert_eq2

    FUNCTION assert_eq3(n1,n2,n3,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        INTEGER, INTENT(IN)          :: n1,n2,n3
        INTEGER                      :: assert_eq3
        IF (n1 == n2 .and. n2 == n3) THEN
            assert_eq3=n1
        ELSE
            WRITE (*,*) 'error: an assert_eq failed with this tag:',string
            WRITE (*,*) 'push any key to finish'
            READ  (*,*)
            STOP 'program terminated by assert_eq3'
        END IF
    END FUNCTION assert_eq3

    FUNCTION assert_eq4(n1,n2,n3,n4,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        INTEGER, INTENT(IN)          :: n1,n2,n3,n4
        INTEGER                      :: assert_eq4
        IF (n1 == n2 .and. n2 == n3 .and. n3 == n4) THEN
            assert_eq4=n1
        ELSE
            WRITE (*,*) 'error: an assert_eq failed with this tag:',string
            WRITE (*,*) 'push any key to finish'
            READ  (*,*)
            STOP 'program terminated by assert_eq4'
        END IF
    END FUNCTION assert_eq4

    FUNCTION assert_eqn(nn,string)
        CHARACTER(LEN=*), INTENT(IN)      :: string
        INTEGER, DIMENSION(:), INTENT(IN) :: nn
        INTEGER                           :: assert_eqn
        IF (all(nn(2:) == nn(1))) THEN
            assert_eqn=nn(1)
        ELSE
            WRITE (*,*) 'error: an assert_eq failed with this tag:',string
            WRITE (*,*) 'push any key to finish'
            READ  (*,*)
            STOP 'program terminated by assert_eqn'
        END IF
    END FUNCTION assert_eqn
!    !
!    !REALLOCATES
!    !
        FUNCTION reallocate_rv(p,n)
	    REAL(SP), DIMENSION(:), POINTER :: p, reallocate_rv
	    INTEGER(I4B), INTENT(IN)        :: n
	    INTEGER(I4B)                    :: nold,ierr
	    ALLOCATE(reallocate_rv(n),stat=ierr)
	    IF (ierr /= 0) CALL error('reallocate_rv: problem in attempt to allocate memory')
	    IF (.not. associated(p)) RETURN
	    nold=size(p)
	    reallocate_rv(1:min(nold,n))=p(1:min(nold,n))
	    DEALLOCATE(p)
	END FUNCTION reallocate_rv
	
	FUNCTION reallocate_dv(p,n)
	    REAL(DP), DIMENSION(:), POINTER :: p, reallocate_dv
	    INTEGER(I4B), INTENT(IN)        :: n
	    INTEGER(I4B)                    :: nold,ierr
	    ALLOCATE(reallocate_dv(n),stat=ierr)
	    IF (ierr /= 0) CALL error('reallocate_dv: problem in attempt to allocate memory')
	    IF (.not. associated(p)) RETURN
	    nold=size(p)
	    reallocate_dv(1:min(nold,n))=p(1:min(nold,n))
	    DEALLOCATE(p)
	END FUNCTION reallocate_dv
    
	FUNCTION reallocate_iv(p,n)
	    INTEGER(I4B), DIMENSION(:), POINTER :: p, reallocate_iv
	    INTEGER(I4B), INTENT(IN)            :: n
	    INTEGER(I4B)                        :: nold,ierr
	    ALLOCATE(reallocate_iv(n),stat=ierr)
	    IF (ierr /= 0) CALL error('reallocate_iv: problem in attempt to allocate memory')
	    IF (.not. associated(p)) RETURN
	    nold=size(p)
	    reallocate_iv(1:min(nold,n))=p(1:min(nold,n))
	    DEALLOCATE(p)
	END FUNCTION reallocate_iv
    
	FUNCTION reallocate_hv(p,n)
	    CHARACTER(1), DIMENSION(:), POINTER :: p, reallocate_hv
	    INTEGER(I4B), INTENT(IN)            :: n
	    INTEGER(I4B)                        :: nold,ierr
	    ALLOCATE(reallocate_hv(n),stat=ierr)
	    IF (ierr /= 0) CALL error('reallocate_hv: problem in attempt to allocate memory')
	    IF (.not. associated(p)) RETURN
	    nold=size(p)
	    reallocate_hv(1:min(nold,n))=p(1:min(nold,n))
	    DEALLOCATE(p)
	END FUNCTION reallocate_hv
	
	FUNCTION reallocate_rm(p,n,m)
	    REAL(SP), DIMENSION(:,:), POINTER :: p, reallocate_rm
	    INTEGER(I4B), INTENT(IN)          :: n,m
	    INTEGER(I4B)                      :: nold,mold,ierr
	    ALLOCATE(reallocate_rm(n,m),stat=ierr)
	    IF (ierr /= 0) CALL error('reallocate_rm: problem in attempt to allocate memory')
	    IF (.not. associated(p)) RETURN
	    nold=size(p,1)
	    mold=size(p,2)
	    reallocate_rm(1:min(nold,n),1:min(mold,m))=&
		    p(1:min(nold,n),1:min(mold,m))
	    DEALLOCATE(p)
	END FUNCTION reallocate_rm
	
	FUNCTION reallocate_dm(p,n,m)
	    REAL(DP), DIMENSION(:,:), POINTER :: p, reallocate_dm
	    INTEGER(I4B), INTENT(IN)          :: n,m
	    INTEGER(I4B)                      :: nold,mold,ierr
	    ALLOCATE(reallocate_dm(n,m),stat=ierr)
	    IF (ierr /= 0) CALL error('reallocate_dm: problem in attempt to allocate memory')
	    IF (.not. associated(p)) RETURN
	    nold=size(p,1)
	    mold=size(p,2)
	    reallocate_dm(1:min(nold,n),1:min(mold,m))=&
		    p(1:min(nold,n),1:min(mold,m))
	    DEALLOCATE(p)
	END FUNCTION reallocate_dm
    
	FUNCTION reallocate_im(p,n,m)
	    INTEGER(I4B), DIMENSION(:,:), POINTER :: p, reallocate_im
	    INTEGER(I4B), INTENT(IN)              :: n,m
	    INTEGER(I4B)                          :: nold,mold,ierr
	    ALLOCATE(reallocate_im(n,m),stat=ierr)
	    IF (ierr /= 0) CALL error('reallocate_im: problem in attempt to allocate memory')
	    IF (.not. associated(p)) RETURN
	    nold=size(p,1)
	    mold=size(p,2)
	    reallocate_im(1:min(nold,n),1:min(mold,m))=&
		    p(1:min(nold,n),1:min(mold,m))
	    DEALLOCATE(p)
	END FUNCTION reallocate_im
!    !
    !ERROR
    !	
    SUBROUTINE error(string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        INTEGER(I4B)                 :: lngtstr
        CHARACTER(10)                :: strlngt               
        lngtstr=LEN_TRIM(string)+11
        WRITE(strlngt,'(I10)')lngtstr
        WRITE(*,*)
        WRITE (*,'(A'//TRIM(ADJUSTL(strlngt))//')')' (X) Error '//TRIM(ADJUSTL(string))
        WRITE(*,'(A27)')'     Push any key to finish'
        READ(*,*)
        STOP 'Program terminated by error'
    END SUBROUTINE error    
    !
    !HARMONIC MEAN
    !
    FUNCTION armncmn(a,b)
        USE anisotype
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: a,b
        REAL(DP)             :: armncmn
        IF ((a==0).OR.(b==0)) THEN
            armncmn=0.0
        ELSE
            armncmn= 2.0_dp/(1.0_dp/a+1.0_dp/b)
        END IF
        !
    END FUNCTION armncmn
END MODULE util