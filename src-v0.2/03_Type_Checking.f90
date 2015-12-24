!-----------------------------------------------------------------------------------------------------------------------------------
!Copyright (c) 2013, Alvarez-Villa O.D., Alvarez-Rendon J.P. GOTTA Ingenieria SAS, All Right Reserved.
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
MODULE tpvrblchckng
    !!
    INTERFACE chckdmn
        SUBROUTINE chckdmn(dmn,typdmn,sbrtn)
            USE anisotype
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN)        :: sbrtn
            INTEGER(I4B), INTENT(IN)        :: typdmn
            TYPE(domain), INTENT(IN)        :: dmn
        END SUBROUTINE chckdmn
    END INTERFACE chckdmn
END MODULE tpvrblchckng

!!
!***********************************************************************************************************************************
!**************************************** SUBRUTINA PARA COMPROBAR DOMINIOS ********************************************************
!******************************** @authores: Alvarez-Villa O.D. & Rendon-Alvarez J.P. **********************************************
!***********************************************************************************************************************************
!
!DESCRIPCI�N:
!Esta subrutina comprueba si una variable tipo de geometria ha sido construida adecuadamente.
!
!DATOS DE ENTRADA:
!dmn                : Variable anisotype geometria a comprobar.                                                               
!typdmn             : Precision en la que se ha construido la variable domain:                                                {entero}
!                            0. Entera
!                            1. Real simple
!                            2. Real doble
!sbrtn              : Subrutina en la que se comprueba la variable anisotype. Se usa para generar mensajes de error.          {character}
!!
SUBROUTINE chckdmn(dmn,typdmn,sbrtn)
    !
    USE anisotype
    USE util, ONLY: error
    !
    IMPLICIT NONE
    !
    CHARACTER(*), INTENT(IN)        :: sbrtn
    INTEGER(I4B), INTENT(IN)        :: typdmn
    TYPE(domain), INTENT(IN)        :: dmn
    !
    LOGICAL(LGT)    :: flexist
    !
    !Existencia de la variable anisotype
    IF ((dmn%exist).EQV.(.false.)) CALL error('in '//TRIM(ADJUSTL(sbrtn))//': Domain has not been built')
    IF ((typdmn/=1).AND.(typdmn/=2)) CALL error('in '//TRIM(ADJUSTL(sbrtn))//': non allowed Domain coordinate number format')
    IF ((dmn%infile).EQV.(.false.)) THEN
        !Formato de almacenamiento de coordenadas
        IF ((typdmn==1).AND.((ASSOCIATED(dmn%Xs).EQV.(.false.)).OR.(ASSOCIATED(dmn%Ys)&
            & .EQV.(.false.)).OR.(ASSOCIATED(dmn%Zs).EQV.(.false.))))&
            & CALL error('in '//TRIM(ADJUSTL(sbrtn))//': X, Y or Z attribute of Domain has not been built properly')
        IF ((typdmn==2).AND.((ASSOCIATED(dmn%Xd).EQV.(.false.)).OR.(ASSOCIATED(dmn%Yd)&
            & .EQV.(.false.)).OR.(ASSOCIATED(dmn%Zd).EQV.(.false.))))&
            & CALL error('in '//TRIM(ADJUSTL(sbrtn))//': X, Y or Z attribute of Domain has not been built properly')
        !Dimensiones de arreglos de almacenamiento
        IF ((typdmn==1).AND.((SIZE(dmn%Xs,1)/=dmn%clmns+1).OR.(SIZE(dmn%Ys,1)/=dmn%rws+1)&
            & .OR.(SIZE(dmn%Zs,1)/=dmn%lvls+1))) CALL error('in '//TRIM(ADJUSTL(sbrtn))//&
            & ': one Domain coordinate array dimension is not appropiate')
        IF ((typdmn==2).AND.((SIZE(dmn%Xd,1)/=dmn%clmns+1).OR.(SIZE(dmn%Yd,1)/=dmn%rws+1)&
            & .OR.(SIZE(dmn%Zd,1)/=dmn%lvls+1))) CALL error('in '//TRIM(ADJUSTL(sbrtn))//&
            & ': one Domain coordinate array dimension is not appropiate')
    ELSE
        INQUIRE(FILE=TRIM(ADJUSTL(dmn%rtdom)),EXIST=flexist)
        IF (.NOT.flexist) &
            & CALL error('in '//TRIM(ADJUSTL(sbrtn))//': Domain file (*.domnRST) does not exist or is not in the correct folder')
        IF (dmn%flind==-9999) CALL error('in'//TRIM(ADJUSTL(sbrtn))//': there is not Domain unit number')
    END IF
    !
END SUBROUTINE chckdmn
!! 