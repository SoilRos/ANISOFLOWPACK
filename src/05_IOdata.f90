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
MODULE iodata
        
    INTERFACE imprtdttrbn
        SUBROUTINE imprtdttrbn(prjctrtd,prjctrt,forminI,forminR,dmngmtry)
            USE anisotype
            USE util
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN)      :: prjctrtd,prjctrt,forminI,forminR    
            TYPE(domain), INTENT(IN)      :: dmngmtry
        END SUBROUTINE imprtdttrbn
    END INTERFACE imprtdttrbn
    
    INTERFACE imprtdtstbn
        SUBROUTINE imprtdtstbn(prjctrtd,prjctrt,forminI,forminR,dmngmtry)
            USE anisotype
            USE util
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN)      :: prjctrtd,prjctrt,forminI,forminR    
            TYPE(domain), INTENT(IN)      :: dmngmtry
        END SUBROUTINE imprtdtstbn
    END INTERFACE imprtdtstbn
    
    INTERFACE aniso_imprtdttrbn
        SUBROUTINE aniso_imprtdttrbn(prjctrtd,prjctrt,forminI,forminR,dmngmtry)
            USE anisotype
            USE util
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN)      :: prjctrtd,prjctrt,forminI,forminR    
            TYPE(domain), INTENT(IN)      :: dmngmtry
        END SUBROUTINE aniso_imprtdttrbn
    END INTERFACE aniso_imprtdttrbn
        
    INTERFACE aniso_imprtdtstbn
        SUBROUTINE aniso_imprtdtstbn(prjctrtd,prjctrt,forminI,forminR,dmngmtry)
            USE anisotype
            USE util
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN)      :: prjctrtd,prjctrt,forminI,forminR    
            TYPE(domain), INTENT(IN)      :: dmngmtry
        END SUBROUTINE aniso_imprtdtstbn
    END INTERFACE aniso_imprtdtstbn
!        
    INTERFACE gtvlr_bn
        SUBROUTINE gtvlr_intbn(u,dmngmtry,i,j,k,vl)
            USE anisotype
            IMPLICIT NONE        
            INTEGER(I4B), INTENT(IN)  :: u,i,j,k
            TYPE(domain), INTENT(IN)  :: dmngmtry
            INTEGER(I4B), INTENT(OUT) :: vl
	END SUBROUTINE gtvlr_intbn
        
	SUBROUTINE gtvlr_dpbn(u,dmngmtry,i,j,k,vl)
            USE anisotype
            IMPLICIT NONE
            INTEGER(I4B), INTENT(IN) :: u,i,j,k
            TYPE(domain), INTENT(IN) :: dmngmtry
            REAL(DP), INTENT(OUT)    :: vl
	END SUBROUTINE gtvlr_dpbn
    END INTERFACE gtvlr_bn
END MODULE iodata

!DUDAS SOBRE PARAMETROS:
!1. En que archivo se localizan las alturas de nivel impuesto y externo?. Estan incluidas en el archivo de niveles inciales? 
!NOTA: Recordar que niveles iniciales asigna valores a las celdas de alturas impuestas y externas.
!2. Por que el problema en regimen permanente no niveles externos?

!***********************************************************************************************************************************
!**********************PREPROCESAMIENTO DE PARAMETROS PARA REGIMEN TRANSITORIO******************************************************
!***************************@author Álvarez-Villa O.D. & Perez K.*******************************************************************
!***********************************************************************************************************************************    
!DESCRIPCION:
!NOMBRE: Imprtdttr
!Subrutina que transcribe los arcvhivos de parametros a un formato apropiado para el realziar la lectura directa de los valores. 
!Dicha transcripcion se puede hacer a un formato de texto o binario.

!DATOS DE ENTRADA:  
!prjctrtd                 :Ruta completa de la carpeta de datos no formateados.
!prjctrt                  :Ruta completa de la carpeta del proyecto.
!forminI                  :Formato de número de lectura
!forminR                  :Formato de número de escritura
!dmngmtry                 :Variable anisotype de la geometría del acuífero

!COMENTARIOS IMPORTANTES SOBRE PARAMETROS
!'.anpck.acte'! (1)celdas activas, (0)inactivas y (2)niveles impuesto. 
!'.anpck.dphe'! (3)celdas de nivel externo.
!'.anpck.qcte'! Tiene tantos identificadores como zonas de nivel impuesto hayan. Tiene valores para cada celda de nivel impuesto.
!'.anpck.qexe'! Tiene tantos identificadores como zonas de nivel externo hayan. Tiene valores para cada celda de nivel externo.
!'.anpck.rcge'! Tiene tantos identificadores como zonas de recaga hayan.
!'.anpck.wele'! Tiene tantos identificadores como pozos hayan.
!'.anpck.cvte'! Conductividades de bloque (m/d). Tiene valores de conductividades para cada celda activa.
!'.anpck.stge'! Coeficientes de almacanamiento. Tiene valores de coeficiente de almacenamiento para cada celda activa.
!'.anpck.rkye'! Conductividad del lecho de río (m/d). Tiene valores de conductancia para cada celda de nivel externo.  
!'.anpck.rvle'! Longitud de la celda de río (m). Tiene valores para cada celda de nivel externo.
!'.anpck.rvbe'! Ancho de la celda de río (m). Tiene valores para cada celda de nivel externo.
!'.anpck.rdhe'! Espesor de lecho del río (m). Tiene valores para cada celda de nivel externo.
!'.anpck.inhe'! Niveles iniciales en cada bloque (m). Tiene valores para cada celda activa.
!
SUBROUTINE imprtdttrbn(prjctrtd,prjctrt,forminI,forminR,dmngmtry)
    USE anisotype
    !
    IMPLICIT NONE
    !
    CHARACTER(*), INTENT(IN)      :: prjctrtd,prjctrt,forminI,forminR    
    TYPE(domain), INTENT(IN)      :: dmngmtry
    !
    CHARACTER(11), DIMENSION(13)  :: extin
    CHARACTER(10), DIMENSION(13)  :: extout 
    CHARACTER(200)                :: rtin,rtout
    INTEGER(I4B)                  :: i,j,u,v,vali,n
    REAL(DP)                      :: val    
    !
    PARAMETER(u=99,v=71)  
    !
    !Un vector que contiene las extensiones de entrada    
    extin(1)  = 'anpck.acte'!Archivo de celdas activas, inactivas y niveles impuesto.
    extin(2)  = 'anpck.dphe'!Archivo de celdas de nivel externo
    extin(3)  = 'anpck.qcte'!Archivo de zonas de nivel impuesto.
    extin(4)  = 'anpck.qexe'!Archivo de zonas de nivel externo.
    extin(5)  = 'anpck.rcge'!Archivo de zonas homogéneas de recarga.
    extin(6)  = 'anpck.wele'!Archivo de localización de pozos.
    extin(7)  = 'anpck.cvte'!Archivo de conductividades de bloque (m/d).
    extin(8)  = 'anpck.stge'!Archivo de coeficientes de almacanamiento.
    extin(9)  = 'anpck.rkye'!Archivo de conductividad del lecho de río (m/d). Un vector de la misma 
    extin(10) = 'anpck.rvle'!Archivo de longitud de la celda de río (m).
    extin(11) = 'anpck.rvbe'!Archivo de ancho de la celda de río (m).
    extin(12) = 'anpck.rdhe'!Archivo de espesor de lecho del río (m).
    extin(13) = 'anpck.inhe'!Archivo de niveles iniciales en cada bloque (m).              
    !
    !Un vector que contiene las extensiones de salida    
    extout(1)  = 'anpck.act'!Archivo de celdas activas, inactivas y niveles impuesto.
    extout(2)  = 'anpck.dph'!Archivo de celdas de nivel externo.
    extout(3)  = 'anpck.qct'!Archivo de zonas de nivel impuesto.
    extout(4)  = 'anpck.qex'!Archivo de zonas de nivel externo.
    extout(5)  = 'anpck.rcg'!Archivo de zonas homogéneas de recarga.
    extout(6)  = 'anpck.wel'!Archivo de localuización de pozos.
    extout(7)  = 'anpck.cvt'!Archivo de conductividades de bloque (m^2/d).
    extout(8)  = 'anpck.stg'!Archivo de coeficientes de almacanamiento.
    extout(9)  = 'anpck.rky'!Archivo de conductividad del lecho de río (m/d).
    extout(10) = 'anpck.rvl'!Archivo de longitud de la celda de río (m).
    extout(11) = 'anpck.rvb'!Archivo de ancho de la celda de río (m).
    extout(12) = 'anpck.rdh'!Archivo de espesor de lecho del río (m).
    extout(13) = 'anpck.inh'!Archivo de niveles iniciales en cada bloque (m).
    !
    !PENDIENTE: Anexar subrutina que verifique creacion de las variables anisotype necesarias.
    WRITE (*,*) '*---|Procesando el grupo de parametros de estado transitorio en binario|---*'
    n=dmngmtry%lvls*dmngmtry%rws*dmngmtry%clmns
    !Los archivos que tienen identificadores enteros
    DO i=1,6
        rtin=trim(prjctrtd)//'/'//extin(i)
        rtout=trim(prjctrt)//'/'//extout(i)
        OPEN(u, FILE=rtin, STATUS='OLD',ACTION='READ')
        OPEN(v, FILE=rtout, STATUS='REPLACE',ACCESS='DIRECT', FORM='UNFORMATTED', RECL=4*6)
        DO j=1,n     
            READ(u, '(('//forminI//'))')vali
            WRITE(v, rec=j)vali        
        END DO
        CLOSE(u); CLOSE(v)
    END DO
    !
    !Los archivos que contienen reales (Valores de parametros)
    DO i=7,13
        rtin=trim(prjctrtd)//'/'//extin(i)
        rtout=trim(prjctrt)//'/'//extout(i)
        OPEN(u, FILE=rtin, STATUS='OLD',ACTION='READ')
        OPEN(v, FILE=rtout, STATUS='REPLACE',ACCESS='DIRECT', FORM='UNFORMATTED', RECL=8*18)
        DO j=1,n       
            READ(u, '(('//forminR//'))')val
            WRITE(v, rec=j)val     
        END DO
        CLOSE(u); CLOSE(v)
    END DO
    WRITE (*,*) '*---|Termina de procesar informacion de estado transitorio en binario|---*'
!    !
END SUBROUTINE imprtdttrbn
!
!Subrutina para preprocesamiento de parametros en estado permanente
!prjctrtd       :Ruta completa de la carpeta de datos no formateados.
!prjctrt        :Ruta completa de la carpeta del proyecto.
!prjnm          :Nombre del proyecto de trabajo
!forminI        :Formato de número de lectura
!forminR        :Formato de número de escritura
!dmngmtry          :Variable anisotype de la geometría del acuífero
!
SUBROUTINE imprtdtstbn(prjctrtd,prjctrt,forminI,forminR,dmngmtry)
    USE anisotype
    !
    IMPLICIT NONE
    !
    CHARACTER(*), INTENT(IN)      :: prjctrtd,prjctrt,forminI,forminR    
    TYPE(domain), INTENT(IN)      :: dmngmtry
    !
    CHARACTER(12), DIMENSION(11)   :: extin
    CHARACTER(11), DIMENSION(11)   :: extout 
    CHARACTER(200)                :: rtin,rtout
    INTEGER(I4B)                  :: i,j,u,v,vali,n
    REAL(DP)                      :: val    
    !
    PARAMETER(u=99,v=60) 
    !
    !Un vector que contiene las extensiones de entrada    
    extin(1) = 'sanpck.acte'!Archivo de celdas activas, inactivas y niveles impuesto permanente.
    extin(2) = 'sanpck.dphe'!Archivo de celdas de nivel externo permanente.
    extin(3) = 'sanpck.qcte'!Archivo de zonas de nivel impuesto permanente.
    extin(4) = 'sanpck.qexe'!Archivo de zonas de nivel externo permanente.
    extin(5) = 'sanpck.rcge'!Archivo de zonas homogéneas de recarga permanente.
    extin(6) = 'sanpck.wele'!Archivo de localización de pozos permanente.
    extin(7) = 'sanpck.cvte'!Archivo de conductividades de bloque (m/d).
    extin(8) = 'sanpck.rkye'!Archivo de conductividad del lecho de río (m/d).
    extin(9) = 'sanpck.rvle'!Archivo de longitud de la celda de río (m).
    extin(10) = 'sanpck.rvbe'!Archivo de ancho de la celda de río (m).
    extin(11) = 'sanpck.rdhe'!Archivo de espesor de lecho del río (m).
    !
    !Un vector que contiene las extensiones de salida    
    extout(1) = 'sanpck.act'!Archivo de celdas activas, inactivas y niveles impuesto permanente.
    extout(2) = 'sanpck.dph'!Archivo de celdas de nivel externo permanente.
    extout(3) = 'sanpck.qct'!Archivo de zonas de nivel impuesto permanente.
    extout(4) = 'sanpck.qex'!Archivo de zonas de nivel externo permanente.
    extout(5) = 'sanpck.rcg'!Archivo de zonas homogéneas de recarga permanente.
    extout(6) = 'sanpck.wel'!Archivo de localización de pozos permanente.
    extout(7) = 'sanpck.cvt'!Archivo de conductividades de bloque (m/d).
    extout(8) = 'sanpck.rky'!Archivo de conductividad del lecho de río (m/d).
    extout(9) = 'sanpck.rvl'!Archivo de longitud de la celda de río (m).
    extout(10) = 'sanpck.rvb'!Archivo de ancho de la celda de río (m).
    extout(11) = 'sanpck.rdh'!Archivo de espesor de lecho del río (m).
    !
    WRITE (*,*) '*---|Procesando el grupo de parametros de estado permanente en binario|---*'
    n=dmngmtry%lvls*dmngmtry%rws*dmngmtry%clmns
    !Los archivos que tienen identificadores enteros
    DO i=1,6
        rtin=trim(prjctrtd)//'/'//extin(i)
        rtout=trim(prjctrt)//'/'//extout(i)
        OPEN(u, FILE=rtin, STATUS='OLD',ACTION='READ')
        OPEN(v, FILE=rtout, STATUS='REPLACE',ACCESS='DIRECT', FORM='UNFORMATTED', RECL=4*6)
        DO j=1,n     
            READ(u, '(('//forminI//'))')vali
            WRITE(v, rec=j)vali     
        END DO
        CLOSE(u); CLOSE(v)
    END DO
    !
    !Los archivos que contienen reales (Valores de parametros)
    DO i=7,11
        rtin=trim(prjctrtd)//'/'//extin(i)
        rtout=trim(prjctrt)//'/'//extout(i)
        OPEN(u, FILE=rtin, STATUS='OLD',ACTION='READ')
        OPEN(v, FILE=rtout, STATUS='REPLACE',ACCESS='DIRECT', FORM='UNFORMATTED', RECL=8*18)
        DO j=1,n       
            READ(u, '(('//forminR//'))')val
            WRITE(v, rec=j)val 
        END DO 
        CLOSE(u); CLOSE(v)
    END DO
    WRITE (*,*) '*---|Termina de procesar informacion de estado permanente en binario|---*'
!    
END SUBROUTINE imprtdtstbn
!
SUBROUTINE aniso_imprtdttrbn(prjctrtd,prjctrt,forminI,forminR,dmngmtry)
    USE anisotype
    !
    IMPLICIT NONE
    !
    CHARACTER(*), INTENT(IN)      :: prjctrtd,prjctrt,forminI,forminR    
    TYPE(domain), INTENT(IN)      :: dmngmtry
    !
    CHARACTER(11), DIMENSION(16)  :: extin
    CHARACTER(10), DIMENSION(16)  :: extout 
    CHARACTER(200)                :: rtin,rtout
    INTEGER(I4B)                  :: i,j,u,v,vali,n
    REAL(DP)                      :: val    
    !
    PARAMETER(u=99,v=71)  
    !
    !Un vector que contiene las extensiones de entrada    
    extin(1)  = 'anpck.acte'!Archivo de celdas activas, inactivas y niveles impuesto.
    extin(2)  = 'anpck.dphe'!Archivo de celdas de nivel externo
    extin(3)  = 'anpck.qcte'!Archivo de zonas de nivel impuesto.
    extin(4)  = 'anpck.qexe'!Archivo de zonas de nivel externo.
    extin(5)  = 'anpck.rcge'!Archivo de zonas homogéneas de recarga.
    extin(6)  = 'anpck.wele'!Archivo de localización de pozos.
    extin(7)  = 'anpck.alfe'!Archivo de angulo de anisotropia de eje x (Angulo entero entre 0 y 360).
    extin(8)  = 'anpck.bete'!Archivo de angulo de anisotropia de eje y (Angulo entero entre 0 y 360).
    extin(9)  = 'anpck.tete'!Archivo de angulo de anisotropia de eje z (Angulo entero entre 0 y 360).
    extin(10) = 'anpck.cvte'!Archivo de conductividades de bloque (m/d).
    extin(11) = 'anpck.stge'!Archivo de coeficientes de almacanamiento.
    extin(12) = 'anpck.rkye'!Archivo de conductividad del lecho de río (m/d). Un vector de la misma 
    extin(13) = 'anpck.rvle'!Archivo de longitud de la celda de río (m).
    extin(14) = 'anpck.rvbe'!Archivo de ancho de la celda de río (m).
    extin(15) = 'anpck.rdhe'!Archivo de espesor de lecho del río (m).
    extin(16) = 'anpck.inhe'!Archivo de niveles iniciales en cada bloque (m).              
    !
    !Un vector que contiene las extensiones de salida    
    extout(1)  = 'anpck.act'!Archivo de celdas activas, inactivas y niveles impuesto.
    extout(2)  = 'anpck.dph'!Archivo de celdas de nivel externo.
    extout(3)  = 'anpck.qct'!Archivo de zonas de nivel impuesto.
    extout(4)  = 'anpck.qex'!Archivo de zonas de nivel externo.
    extout(5)  = 'anpck.rcg'!Archivo de zonas homogéneas de recarga.
    extout(6)  = 'anpck.wel'!Archivo de localuización de pozos.
    extout(7)  = 'anpck.alf'!Archivo de angulo de anisotropia de eje x (Angulo entero entre 0 y 360).
    extout(8)  = 'anpck.bet'!Archivo de angulo de anisotropia de eje y (Angulo entero entre 0 y 360).
    extout(9)  = 'anpck.tet'!Archivo de angulo de anisotropia de eje z (Angulo entero entre 0 y 360).
    extout(10) = 'anpck.cvt'!Archivo de conductividades de bloque (m^2/d).
    extout(11) = 'anpck.stg'!Archivo de coeficientes de almacanamiento.
    extout(12) = 'anpck.rky'!Archivo de conductividad del lecho de río (m/d).
    extout(13) = 'anpck.rvl'!Archivo de longitud de la celda de río (m).
    extout(14) = 'anpck.rvb'!Archivo de ancho de la celda de río (m).
    extout(15) = 'anpck.rdh'!Archivo de espesor de lecho del río (m).
    extout(16) = 'anpck.inh'!Archivo de niveles iniciales en cada bloque (m).
    !
    !PENDIENTE: Anexar subrutina que verifique creacion de las variables anisotype necesarias.
    WRITE (*,*) '*---|Procesando el grupo de parametros de estado transitorio en binario|---*'
    n=dmngmtry%lvls*dmngmtry%rws*dmngmtry%clmns
    !Los archivos que tienen identificadores enteros
    DO i=1,9
        rtin=trim(prjctrtd)//'/'//extin(i)
        rtout=trim(prjctrt)//'/'//extout(i)
        OPEN(u, FILE=rtin, STATUS='OLD',ACTION='READ')
        OPEN(v, FILE=rtout, STATUS='REPLACE',ACCESS='DIRECT', FORM='UNFORMATTED', RECL=4*6)
        DO j=1,n     
            READ(u, '(('//forminI//'))')vali
            WRITE(v, rec=j)vali        
        END DO
        CLOSE(u); CLOSE(v)
    END DO
    !
    !Los archivos que contienen reales (Valores de parametros)
    DO i=10,16
        rtin=trim(prjctrtd)//'/'//extin(i)
        rtout=trim(prjctrt)//'/'//extout(i)
        OPEN(u, FILE=rtin, STATUS='OLD',ACTION='READ')
        OPEN(v, FILE=rtout, STATUS='REPLACE',ACCESS='DIRECT', FORM='UNFORMATTED', RECL=8*18)
        DO j=1,n       
            READ(u, '(('//forminR//'))')val
            WRITE(v, rec=j)val     
        END DO
        CLOSE(u); CLOSE(v)
    END DO
    WRITE (*,*) '*---|Termina de procesar informacion de estado transitorio en binario|---*'
!    !
END SUBROUTINE aniso_imprtdttrbn
!
SUBROUTINE aniso_imprtdtstbn(prjctrtd,prjctrt,forminI,forminR,dmngmtry)
    USE anisotype
    !
    IMPLICIT NONE
    !
    CHARACTER(*), INTENT(IN)      :: prjctrtd,prjctrt,forminI,forminR    
    TYPE(domain), INTENT(IN)      :: dmngmtry
    !
    CHARACTER(12), DIMENSION(14)   :: extin
    CHARACTER(11), DIMENSION(14)   :: extout 
    CHARACTER(200)                :: rtin,rtout
    INTEGER(I4B)                  :: i,j,u,v,vali,n
    REAL(DP)                      :: val    
    !
    PARAMETER(u=99,v=60) 
    !
    !Un vector que contiene las extensiones de entrada    
    extin(1) = 'sanpck.acte'!Archivo de celdas activas, inactivas y niveles impuesto permanente.
    extin(2) = 'sanpck.dphe'!Archivo de celdas de nivel externo permanente.
    extin(3) = 'sanpck.qcte'!Archivo de zonas de nivel impuesto permanente.
    extin(4) = 'sanpck.qexe'!Archivo de zonas de nivel externo permanente.
    extin(5) = 'sanpck.rcge'!Archivo de zonas homogéneas de recarga permanente.
    extin(6) = 'sanpck.wele'!Archivo de localización de pozos permanente.
    extin(7) = 'sanpck.alfe'!Archivo de angulo de anisotropia de eje x (Angulo entero entre 0 y 360).
    extin(8) = 'sanpck.bete'!Archivo de angulo de anisotropia de eje y (Angulo entero entre 0 y 360).
    extin(9) = 'sanpck.tete'!Archivo de angulo de anisotropia de eje z (Angulo entero entre 0 y 360).
    extin(10) = 'sanpck.cvte'!Archivo de conductividades de bloque (m/d).
    extin(11) = 'sanpck.rkye'!Archivo de conductividad del lecho de río (m/d).
    extin(12) = 'sanpck.rvle'!Archivo de longitud de la celda de río (m).
    extin(13) = 'sanpck.rvbe'!Archivo de ancho de la celda de río (m).
    extin(14) = 'sanpck.rdhe'!Archivo de espesor de lecho del río (m).
    !
    !Un vector que contiene las extensiones de salida    
    extout(1) = 'sanpck.act'!Archivo de celdas activas, inactivas y niveles impuesto permanente.
    extout(2) = 'sanpck.dph'!Archivo de celdas de nivel externo permanente.
    extout(3) = 'sanpck.qct'!Archivo de zonas de nivel impuesto permanente.
    extout(4) = 'sanpck.qex'!Archivo de zonas de nivel externo permanente.
    extout(5) = 'sanpck.rcg'!Archivo de zonas homogéneas de recarga permanente.
    extout(6) = 'sanpck.wel'!Archivo de localización de pozos permanente.
    extout(7) = 'sanpck.alf'!Archivo de angulos de anisotropia de eje x (Angulo entero entre 0 y 360).
    extout(8) = 'sanpck.bet'!Archivo de angulos de anisotropia de eje y (Angulo entero entre 0 y 360).
    extout(9) = 'sanpck.tet'!Archivo de angulos de anisotropia de eje z (Angulo entero entre 0 y 360).
    extout(10) = 'sanpck.cvt'!Archivo de conductividades de bloque (m/d).
    extout(11) = 'sanpck.rky'!Archivo de conductividad del lecho de río (m/d).
    extout(12) = 'sanpck.rvl'!Archivo de longitud de la celda de río (m).
    extout(13) = 'sanpck.rvb'!Archivo de ancho de la celda de río (m).
    extout(14) = 'sanpck.rdh'!Archivo de espesor de lecho del río (m).
    !
    WRITE (*,*) '*---|Procesando el grupo de parametros de estado permanente en binario|---*'
    n=dmngmtry%lvls*dmngmtry%rws*dmngmtry%clmns
    !Los archivos que tienen identificadores enteros
    DO i=1,9
        rtin=trim(prjctrtd)//'/'//extin(i)
        rtout=trim(prjctrt)//'/'//extout(i)
        OPEN(u, FILE=rtin, STATUS='OLD',ACTION='READ')
        OPEN(v, FILE=rtout, STATUS='REPLACE',ACCESS='DIRECT', FORM='UNFORMATTED', RECL=4*6)
        DO j=1,n     
            READ(u, '(('//forminI//'))')vali
            WRITE(v, rec=j)vali             
        END DO
        CLOSE(u); CLOSE(v)
    END DO
    !
    !Los archivos que contienen reales (Valores de parametros)
    DO i=10,14
        rtin=trim(prjctrtd)//'/'//extin(i)
        rtout=trim(prjctrt)//'/'//extout(i)
        OPEN(u, FILE=rtin, STATUS='OLD',ACTION='READ')
        OPEN(v, FILE=rtout, STATUS='REPLACE',ACCESS='DIRECT', FORM='UNFORMATTED', RECL=8*18)
        DO j=1,n       
            READ(u, '(('//forminR//'))')val
            WRITE(v, rec=j)val 
        END DO 
        CLOSE(u); CLOSE(v)
    END DO
    WRITE (*,*) '*---|Termina de procesar informacion de estado permanente en binario|---*'
!    
END SUBROUTINE aniso_imprtdtstbn
!!
!***********************************************************************************************************************************
!**********************GET VALORES**************************************************************************************************
!****************@author Álvarez-Villa O.D. & Perez K.******************************************************************************
!***********************************************************************************************************************************
!
!Una rutina que lee el valor de una vector que representa una matriz
!tridimensional en un archivo de parámetros (Similar al método getValor, HSJ). 
!Contiene subrutinas para lectura de un archivo binario.
!Los parámetros que recibe son los siguientes:
!u          :Posicion en memoria del archivo de datos
!i          :Fila del raster
!j          :Columna del raster
!k          :Capa del raster
!vl         :Valor a devolver
!
!Grupo de funciones get para lectura de archivos en binario
!
SUBROUTINE gtvlr_intbn(u,dmngmtry,i,j,k,vl)
    USE anisotype
    USE util, ONLY : error
    !
    IMPLICIT NONE        
    !
    INTEGER(I4B), INTENT(IN)  :: u,i,j,k
    TYPE(domain), INTENT(IN)  :: dmngmtry
    INTEGER(I4B), INTENT(OUT) :: vl
    !
    INTEGER(I4B) :: pstn=0    
    CHARACTER(1) :: intro=' '
    !
    pstn=0; vl=0
    pstn=dmngmtry%clmns*dmngmtry%rws*(k-1)+dmngmtry%clmns*(j-1)+i
    !
    !Checkeo
    IF(pstn>dmngmtry%rws*dmngmtry%clmns*dmngmtry%lvls)CALL error('Posicion solicitada esta fuera del dominio espacial: gtvlr_int ')
    !Valor obtenido
    READ(u,REC=pstn)vl
    !
END SUBROUTINE gtvlr_intbn
!
SUBROUTINE gtvlr_dpbn(u,dmngmtry,i,j,k,vl)
    USE anisotype
    USE util, ONLY : error
    !
    IMPLICIT NONE
    !
    INTEGER(I4B), INTENT(IN)  :: u,i,j,k
    TYPE(domain), INTENT(IN)  :: dmngmtry
    REAL(DP), INTENT(OUT)     :: vl
    !
    INTEGER(I4B) :: pstn=0  
    !
    pstn=dmngmtry%clmns*dmngmtry%rws*(k-1)+dmngmtry%clmns*(j-1)+i
    !
!    !Checkeo
    IF(pstn>dmngmtry%rws*dmngmtry%clmns*dmngmtry%lvls)CALL error('Posicion solicitada esta fuera d dominio espacial:gtvlr_dp')
    !Valor obtenido
    READ(u,REC=pstn)vl    
    !
END SUBROUTINE gtvlr_dpbn