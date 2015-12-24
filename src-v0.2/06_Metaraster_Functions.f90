!-----------------------------------------------------------------------------------------------------------------------------------
!Copyright (c) 2013, Alvarez-Villa O.D. & Rend�n-�lvarez & Universidad Politecnica de Valencia & GOTTA Ingenieria SAS, All Right Reserved.
!
!This Program is distributed in the hope that they will be used but WITHOUT ANY WARRANTY. No autor or distributor accepts 
!responsability to anyone for consequences of using them or for whether they serve !any particula porpouse or work at all, unless he
!says so in writing. Everyone is granted permission to copy, modify an distribute this program, but only under the conditions that 
!this notice and above copyright notice remains intact.
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!Este modulo contiene rutinas para manipular los objetos raster 3D. Las rutinas soportadas son:
!1. Constructor y destructor de la variable tipo metaraster(type) desde archivo
!2. Constructor y destructor de la variable tipo metaraster directamente
!3. Rutina para manipular los metadatos del raster
!4. Rutinas de entradas y salida para escribir y leer los raster
!-----------------------------------------------------------------------------------------------------------------------------------
!
MODULE raster    
    !!
    INTERFACE bld_mtrstrff
        SUBROUTINE bld_mtrstrff_i(prjctrt,flnm,vrbl,units,X,Y,Z,dataval,mtrstr)
            USE festellustype
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN)                                        :: prjctrt,flnm
            CHARACTER(10), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(OUT)   :: units
            CHARACTER(20), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(OUT)   :: vrbl
            INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE, TARGET, INTENT(OUT)  :: dataval
            REAL(SP), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(OUT)        :: X,Y,Z
            TYPE(metaraster), INTENT(OUT)                                   :: mtrstr
        END SUBROUTINE bld_mtrstrff_i
        !BL
        SUBROUTINE bld_mtrstrff_sp(prjctrt,flnm,vrbl,units,X,Y,Z,dataval,mtrstr)
            USE festellustype
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN)                                        :: prjctrt,flnm
            CHARACTER(10), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(OUT)   :: units
            CHARACTER(20), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(OUT)   :: vrbl
            REAL(SP), DIMENSION(:,:), ALLOCATABLE, TARGET, INTENT(OUT)      :: dataval
            REAL(SP), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(OUT)        :: X,Y,Z
            TYPE(metaraster), INTENT(OUT)                                   :: mtrstr
        END SUBROUTINE bld_mtrstrff_sp
        !BL
        SUBROUTINE bld_mtrstrff_dp(prjctrt,flnm,vrbl,units,X,Y,Z,dataval,mtrstr)
            USE festellustype
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN)                                        :: prjctrt,flnm
            CHARACTER(10), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(OUT)   :: units
            CHARACTER(20), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(OUT)   :: vrbl
            REAL(DP), DIMENSION(:,:), ALLOCATABLE, TARGET, INTENT(OUT)      :: dataval
            REAL(DP), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(OUT)        :: X,Y,Z
            TYPE(metaraster), INTENT(OUT)                                   :: mtrstr
        END SUBROUTINE bld_mtrstrff_dp
    END INTERFACE bld_mtrstrff
    !!
    INTERFACE bld_mtrstrif
        SUBROUTINE bld_mtrstrif(prjctrt,flnm,typdmn,typnum,domflind,datflind,mtrstr)
            USE festellustype
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN)      :: prjctrt,flnm
            INTEGER(I4B), INTENT(IN)      :: typdmn,typnum,domflind,datflind
            TYPE(metaraster), INTENT(OUT) :: mtrstr
        END SUBROUTINE bld_mtrstrif
    END INTERFACE bld_mtrstrif
    !!
    INTERFACE bld_mtrstr
        SUBROUTINE bld_mtrstr_i(falt,name,nvbls,vrbl,units,rgn,info,domn,dataval,mtrstr)    
            USE festellustype
            IMPLICIT NONE
            CHARACTER(20), INTENT(IN)                                       :: name,rgn
            CHARACTER(10), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(IN)    :: units
            CHARACTER(20), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(IN)    :: vrbl
            CHARACTER(500), INTENT(IN)                                      :: info
            INTEGER(I4B), INTENT(IN)                                        :: falt,nvbls
            INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE, TARGET, INTENT(IN)   :: dataval
            TYPE(domain), INTENT(IN)                                        :: domn
            TYPE(metaraster), INTENT(OUT)                                   :: mtrstr
        END SUBROUTINE bld_mtrstr_i
        !BL
        SUBROUTINE bld_mtrstr_sp(falt,name,nvbls,vrbl,units,rgn,info,domn,dataval,mtrstr)    
            USE festellustype
            IMPLICIT NONE
            CHARACTER(20), INTENT(IN)                                       :: name,rgn
            CHARACTER(10), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(IN)    :: units
            CHARACTER(20), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(IN)    :: vrbl
            CHARACTER(500), INTENT(IN)                                      :: info
            INTEGER(I4B), INTENT(IN)                                        :: nvbls
            REAL(SP), INTENT(IN)                                            :: falt
            REAL(SP), DIMENSION(:,:), ALLOCATABLE, TARGET, INTENT(IN)       :: dataval
            TYPE(domain), INTENT(IN)                                        :: domn
            TYPE(metaraster), INTENT(OUT)                                   :: mtrstr
        END SUBROUTINE bld_mtrstr_sp
        !BL
        SUBROUTINE bld_mtrstr_dp(falt,name,nvbls,vrbl,units,rgn,info,domn,dataval,mtrstr)    
            USE festellustype
            IMPLICIT NONE
            CHARACTER(20), INTENT(IN)                                       :: name,rgn
            CHARACTER(10), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(IN)    :: units
            CHARACTER(20), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(IN)    :: vrbl
            CHARACTER(500), INTENT(IN)                                      :: info
            INTEGER(I4B), INTENT(IN)                                        :: nvbls
            REAL(DP), INTENT(IN)                                            :: falt
            REAL(DP), DIMENSION(:,:), ALLOCATABLE, TARGET, INTENT(IN)       :: dataval
            TYPE(domain), INTENT(IN)                                        :: domn
            TYPE(metaraster), INTENT(OUT)                                   :: mtrstr
        END SUBROUTINE bld_mtrstr_dp
    END INTERFACE bld_mtrstr
    !!
    INTERFACE svmetaraster
        SUBROUTINE svmetaraster(prjctrt,fldr,flnm,mtrstr)
            USE festellustype
            IMPLICIT NONE
            CHARACTER(*), INTENT(IN)         :: prjctrt,fldr
            CHARACTER(*), INTENT(INOUT)      :: flnm
            TYPE(metaraster), INTENT(INOUT)  :: mtrstr
        END SUBROUTINE svmetaraster
    END INTERFACE svmetaraster
    !!
    INTERFACE
        SUBROUTINE dstry_metaraster(mtrstr,switch)    
            USE festellustype
            IMPLICIT NONE
            LOGICAL(LGT), INTENT(IN)        :: switch
            TYPE(metaraster), INTENT(INOUT) :: mtrstr
        END SUBROUTINE dstry_metaraster
    END INTERFACE
    !!    
    INTERFACE gt_mtrstrvl
        SUBROUTINE gt_mtrstrvl_i(mtrstr,switch,vl,i,j,k,v,l)
            USE festellustype
            IMPLICIT NONE
            INTEGER(I4B), INTENT(OUT)              :: vl
            INTEGER(I4B), OPTIONAL, INTENT(IN)     :: i,j,k,v,l
            LOGICAL(LGT)                           :: switch
            TYPE(metaraster), INTENT(IN)           :: mtrstr
        END SUBROUTINE gt_mtrstrvl_i
        !BL
        SUBROUTINE gt_mtrstrvl_sp(mtrstr,switch,vl,i,j,k,v,l)
            USE festellustype
            IMPLICIT NONE 
            REAL(SP), INTENT(OUT)                  :: vl
            INTEGER(I4B), OPTIONAL, INTENT(IN)     :: i,j,k,v,l
            LOGICAL(LGT)                           :: switch
            TYPE(metaraster), INTENT(IN)           :: mtrstr
        END SUBROUTINE gt_mtrstrvl_sp
        !BL
        SUBROUTINE gt_mtrstrvl_dp(mtrstr,switch,vl,i,j,k,v,l)
            USE festellustype
            IMPLICIT NONE
            REAL(DP), INTENT(OUT)                  :: vl
            INTEGER(I4B), OPTIONAL, INTENT(IN)     :: i,j,k,v,l
            LOGICAL(LGT)                           :: switch
            TYPE(metaraster), INTENT(IN)           :: mtrstr
        END SUBROUTINE gt_mtrstrvl_dp
    END INTERFACE gt_mtrstrvl
    !!
    INTERFACE smpl_mtrstr_pnts
        SUBROUTINE smpl_mtrstr_pnts_i(mtrstr,vblei,orstval,XYZ,smpldvls)
            USE festellustype
            IMPLICIT NONE
            INTEGER(I4B), INTENT(IN)                         :: vblei,orstval
            INTEGER(I4B), DIMENSION(:), POINTER, INTENT(OUT) :: smpldvls
            REAL(SP), DIMENSION(:,:), INTENT(IN)             :: XYZ
            TYPE(metaraster), INTENT(IN)                     :: mtrstr
        END SUBROUTINE smpl_mtrstr_pnts_i
        !BL
        SUBROUTINE smpl_mtrstr_pnts_sp(mtrstr,vblei,orstval,XYZ,smpldvls)
            USE festellustype
            IMPLICIT NONE
            INTEGER(I4B), INTENT(IN)                         :: vblei
            REAL(SP), INTENT(IN)                             :: orstval
            REAL(SP), DIMENSION(:), POINTER, INTENT(OUT)     :: smpldvls
            REAL(SP), DIMENSION(:,:), INTENT(IN)             :: XYZ
            TYPE(metaraster), INTENT(IN)                     :: mtrstr
        END SUBROUTINE smpl_mtrstr_pnts_sp
        !BL
        SUBROUTINE smpl_mtrstr_pnts_dp(mtrstr,vblei,orstval,XYZ,smpldvls)
            USE festellustype
            IMPLICIT NONE
            INTEGER(I4B), INTENT(IN)                         :: vblei
            REAL(DP), INTENT(IN)                             :: orstval
            REAL(DP), DIMENSION(:), POINTER, INTENT(OUT)     :: smpldvls
            REAL(DP), DIMENSION(:,:), INTENT(IN)             :: XYZ
            TYPE(metaraster), INTENT(IN)                     :: mtrstr
        END SUBROUTINE smpl_mtrstr_pnts_dp
    END INTERFACE smpl_mtrstr_pnts    
    !!
END MODULE raster
!
!***********************************************************************************************************************************
!***************************** IMPORTACI�N DE UNA VARIABLE TIPO METARASTER DESDE ARCHIVO *******************************************
!******************************** @authores �lvarez-Villa O.D. & Rend�n-�lvarez J.P. ***********************************************
!***********************************************************************************************************************************
!
!DESCRIPCI�N:
!Una subrutina que llena los atributos del objeto tipo Metaraster desde tres archivos de texto (*.metaRST, *.domnRST, *.dataRST) 
!almacenados en el disco duro
!
!DATOS DE ENTRADA:
!prjctrt        : Ruta del proyecto de trabajo                                                                            {caracter}
!flnm           : Nombre del proyecto de trabajo                                                                          {caracter}
!
!DATOS DE SALIDA:
!vrbl           : Arreglo para almacenar los nombres de las variables                                               {caracter}(:)[O]
!units          : Arreglo para almacenar las unidades de las variables                                              {caracter}(:)[O]
!X,Y,Z          : Arreglos para almacenar las coordenadas del dominio espacial                                          {real}(:)[O]
!dataval        : Arreglo para almacenar los valores de las variables del Metaraster                         {real � entero}(:,:)[O]
!mtrstr         : Variable tipo Metaraster a la cual se le asignan los par�metros                                       {metaraster}
!
SUBROUTINE bld_mtrstrff_i(prjctrt,flnm,vrbl,units,X,Y,Z,dataval,mtrstr)
    !
    USE festellustype
    USE util,   ONLY : error
    USE domainprop, ONLY : bld_domainff
    USE raster, ONLY :  bld_mtrstr
    !
    IMPLICIT NONE
    !
    CHARACTER(*), INTENT(IN)                                        :: prjctrt,flnm
    CHARACTER(10), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(OUT)   :: units
    CHARACTER(20), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(OUT)   :: vrbl
    INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE, TARGET, INTENT(OUT)  :: dataval
    REAL(SP), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(OUT)        :: X,Y,Z
    TYPE(metaraster), INTENT(OUT)                                   :: mtrstr
    !
    CHARACTER(1)                                        :: enter
    CHARACTER(2)                                        :: rtdlmtng
    CHARACTER(20)                                       :: name,rgn
    CHARACTER(200)                                      :: rtm,rtd
    CHARACTER(500)                                      :: info
    INTEGER(I4B)                                        :: um,ud,i,j,nvbls,falt,hrs,mnts
    REAL(DP)                                            :: timea,timeb,exetime,scnds
    LOGICAL(LGT)                                        :: flexist
    TYPE(domain)                                        :: dmngmtry
    !
    PARAMETER(um=30,ud=31)
    !
    CALL CPU_TIME(timea)                                                                                                            Se pregunta por la hora para controlar el tiempo de ejecuci�n
    WRITE(*,'(A59)')'*---|Importing Metaraster (*.metaRST, *.dataRST) files|---*'                                                   Impresi�n del aviso de inicio de la importaci�n
    !Definici�n del caracter de delimitaci�n de rutas de acuerdo con el sistema operativo
    IF (INDEX(prjctrt,'/')==0) THEN
        rtdlmtng='\\'                                                                                                               Delimitador en MS-DOS
    ELSE
        rtdlmtng='//'                                                                                                               Delimitador en Unix-OS
    END IF
    !Verificaci�n de la existencia y pertinencia de los datos de entrada
    IF (LEN_TRIM(TRIM(ADJUSTL(prjctrt))//rtdlmtng//'Data'//rtdlmtng//'RST'//rtdlmtng//TRIM(ADJUSTL(flnm))//'.metaRST') > 200) THEN
        CALL error('in bld_mtrstrff_i: the path of the Metaraster is very long')    
    ELSE
        rtm=TRIM(ADJUSTL(prjctrt))//rtdlmtng//'Data'//rtdlmtng//'RST'//rtdlmtng//TRIM(ADJUSTL(flnm))//'.metaRST'                    Archivo en el que  se almacenan los metadatos
        rtd=TRIM(ADJUSTL(prjctrt))//rtdlmtng//'Data'//rtdlmtng//'RST'//rtdlmtng//TRIM(ADJUSTL(flnm))//'.dataRST'                    Archivo en el que  se almacena el arreglo de datos
    END IF
    INQUIRE(FILE=TRIM(ADJUSTL(rtm)),EXIST=flexist)
    IF (.NOT.flexist) CALL error('in bld_mtrstrff_i: Metaraster metadata file does not exist')                                      Verifica que el archivo con los metadatos exista y est� en la ruta correcta
    INQUIRE(FILE=TRIM(ADJUSTL(rtd)),EXIST=flexist)
    IF (.NOT.flexist) CALL error('in bld_mtrstrff_i: Metaraster data file does not exist')                                          Verifica que el archivo con el arreglo de datos exista y est� en la ruta correcta
    !Lectura del archivo de metadatos
    OPEN(um,FILE=rtm,STATUS='OLD')                                                                                                  Abre el archivo de los metadatos
    READ(um,*);READ(um,'(A20)')name                                                                                                 Nombre del metaraster                                                                                    Entero que indica el tipo de n�mero para manejar los datos
    READ(um,*);READ(um,*);READ(um,'(I9)')falt                                                                                       Valor de datos faltantes en precisi�n entera
    READ(um,*);READ(um,*);READ(um,'(I9)')nvbls                                                                                      N�mero de variables
    ALLOCATE(vrbl(nvbls),units(nvbls))                                                                                              Aloja los arreglos de variables y unidades
    READ(um,*);READ(um,*)   
    DO i=1,nvbls
        READ(um,'(A20)')vrbl(i)                                                                                                     Nombres de la variables almacenadas en el Metaraster
    END DO
    READ(um,*);READ(um,*)
    DO i=1,nvbls
        READ(um,'(A10)')units(i)                                                                                                    Unidades de las variables almacenadas en el Metaraster
    END DO
    READ(um,*);READ(um,*);READ(um,'(A20)')rgn                                                                                       Regi�n en la que se define el Metaraster
    READ(um,*);READ(um,*);READ(um,'(A500)')info                                                                                     Informaci�n de inter�s sobre el Metaraster
    CLOSE(um)                                                                                                                       Cierra el archivo de los metadatos
    !Construye el dominio espacial
    CALL bld_domainff(prjctrt,flnm,X,Y,Z,dmngmtry)
    !Lectura del archivo que almacena el arreglo de datos
    ALLOCATE(dataval(dmngmtry%clmns*dmngmtry%rws*dmngmtry%lvls,nvbls))                                                              Aloja el arreglo de datos con precisi�n entera
    OPEN(ud,FILE=rtd,STATUS='OLD')                                                                                                  Abre el archivo de datos
    DO j=1,nvbls
        DO i=1,dmngmtry%clmns*dmngmtry%rws*dmngmtry%lvls
            READ(ud,'(SP,I11,A1)')dataval(i,j),enter                                                                                Y lo llena con la informaci�n del archivo de datos
        END DO
    END DO
    CALL bld_mtrstr(falt,name,nvbls,vrbl,units,rgn,info,dmngmtry,dataval,mtrstr)                                                    Finalmente construye el Metaraster
    CLOSE(ud)                                                                                                                       Cierra el archivo de datos
    CALL CPU_TIME(timeb)                                                                                                            Pregunta por la hora para estimar el tiempo de ejecuci�n del programa
    exetime=timeb-timea                                                                                                             Calcula el tiempo de ejecuci�n
    hrs=FLOOR(exetime/3600);mnts=FLOOR((exetime/3600-hrs)*60);scnds=exetime-mnts*60-hrs*3600                                        Convierte las unidades del tiempo de ejecuci�n
    WRITE(*,'(A101,2(I2,A2),F4.1,A6)')'*---|Importing Metaraster (*.metaRST, *.dataRST) files has been finished.'//&                Imprime el mensaje de finalizaci�n satisfactoria del programa
        & '            Execution time: ',hrs,'h ',mnts,'m ',scnds,'s|---*'
    !
END SUBROUTINE bld_mtrstrff_i
!BL
SUBROUTINE bld_mtrstrff_sp(prjctrt,flnm,vrbl,units,X,Y,Z,dataval,mtrstr)
    !
    USE festellustype
    USE util,   ONLY : error
    USE iofrmtng, ONLY : strtonum
    USE domainprop, ONLY : bld_domainff
    USE raster, ONLY :  bld_mtrstr
    !
    IMPLICIT NONE
    !
    CHARACTER(*), INTENT(IN)                                        :: prjctrt,flnm
    CHARACTER(10), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(OUT)   :: units
    CHARACTER(20), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(OUT)   :: vrbl
    REAL(SP), DIMENSION(:,:), ALLOCATABLE, TARGET, INTENT(OUT)      :: dataval
    REAL(SP), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(OUT)        :: X,Y,Z
    TYPE(metaraster), INTENT(OUT)                                   :: mtrstr
    !
    CHARACTER(1)                                        :: enter
    CHARACTER(2)                                        :: rtdlmtng
    CHARACTER(20)                                       :: name,rgn,strnum
    CHARACTER(200)                                      :: rtm,rtd
    CHARACTER(500)                                      :: info
    INTEGER(I4B)                                        :: um,ud,i,j,nvbls,hrs,mnts
    REAL(SP)                                            :: falt
    REAL(DP)                                            :: timea,timeb,exetime,scnds
    LOGICAL(LGT)                                        :: flexist
    TYPE(domain)                                        :: dmngmtry
    !
    PARAMETER(um=30,ud=31)
    !
    CALL CPU_TIME(timea)                                                                                                            Se pregunta por la hora para controlar el tiempo de ejecuci�n
    WRITE(*,'(A59)')'*---|Importing Metaraster (*.metaRST, *.dataRST) files|---*'                                                   Impresi�n del aviso de inicio de la importaci�n
    !Definici�n del caracter de delimitaci�n de rutas de acuerdo con el sistema operativo
    IF (INDEX(prjctrt,'/')==0) THEN
        rtdlmtng='\\'                                                                                                               Delimitador en MS-DOS
    ELSE
        rtdlmtng='//'                                                                                                               Delimitador en Unix-OS
    END IF
    !Verificaci�n de la existencia y pertinencia de los datos de entrada
    IF (LEN_TRIM(TRIM(ADJUSTL(prjctrt))//rtdlmtng//'Data'//rtdlmtng//'RST'//rtdlmtng//TRIM(ADJUSTL(flnm))//'.metaRST') > 200) THEN
        CALL error('in bld_mtrstrff_sp: the path of the Metaraster is very long')    
    ELSE
        rtm=TRIM(ADJUSTL(prjctrt))//rtdlmtng//'Data'//rtdlmtng//'RST'//rtdlmtng//TRIM(ADJUSTL(flnm))//'.metaRST'                    Archivo en el que  se almacenan los metadatos
        rtd=TRIM(ADJUSTL(prjctrt))//rtdlmtng//'Data'//rtdlmtng//'RST'//rtdlmtng//TRIM(ADJUSTL(flnm))//'.dataRST'                    Archivo en el que  se almacena el arreglo de datos
    END IF
    INQUIRE(FILE=TRIM(ADJUSTL(rtm)),EXIST=flexist)
    IF (.NOT.flexist) CALL error('in bld_mtrstrff_sp: Metaraster metadata file does not exist')                                     Verifica que el archivo con los metadatos exista y est� en la ruta correcta
    INQUIRE(FILE=TRIM(ADJUSTL(rtd)),EXIST=flexist)
    IF (.NOT.flexist) CALL error('in bld_mtrstrff_sp: Metaraster data file does not exist')                                         Verifica que el archivo con el arreglo de datos exista y est� en la ruta correcta
    !Lectura del archivo de metadatos
    OPEN(um,FILE=rtm,STATUS='OLD')                                                                                                  Abre el archivo de los metadatos
    READ(um,*);READ(um,'(A20)')name                                                                                                 Nombre del metaraster                                                                                    Entero que indica el tipo de n�mero para manejar los datos
    READ(um,*);READ(um,*);READ(um,'(A20)')strnum                                                                                  
    CALL strtonum(strnum,falt)                                                                                                      Valor de datos faltantes en precisi�n real simple
    READ(um,*);READ(um,*);READ(um,'(I9)')nvbls                                                                                      N�mero de variables
    ALLOCATE(vrbl(nvbls),units(nvbls))                                                                                              Aloja los arreglos de variables y unidades
    READ(um,*);READ(um,*)   
    DO i=1,nvbls
        READ(um,'(A20)')vrbl(i)                                                                                                     Nombres de la variables almacenadas en el Metaraster
    END DO
    READ(um,*);READ(um,*)
    DO i=1,nvbls
        READ(um,'(A10)')units(i)                                                                                                    Unidades de las variables almacenadas en el Metaraster
    END DO
    READ(um,*);READ(um,*);READ(um,'(A20)')rgn                                                                                       Regi�n en la que se define el Metaraster
    READ(um,*);READ(um,*);READ(um,'(A500)')info                                                                                     Informaci�n de inter�s sobre el Metaraster
    CLOSE(um)                                                                                                                       Cierra el archivo de los metadatos
    !Construye el dominio espacial
    CALL bld_domainff(prjctrt,flnm,X,Y,Z,dmngmtry)
    !Lectura del archivo que almacena el arreglo de datos
    ALLOCATE(dataval(dmngmtry%clmns*dmngmtry%rws*dmngmtry%lvls,nvbls))                                                              Aloja el arreglo de datos con precisi�n real simple
    OPEN(ud,FILE=rtd,STATUS='OLD')                                                                                                  Abre el archivo de datos
    DO j=1,nvbls
        DO i=1,dmngmtry%clmns*dmngmtry%rws*dmngmtry%lvls
            READ(ud,'(SP,ES18.11E2,A1)')dataval(i,j),enter                                                                          Y lo llena con la informaci�n del archivo de datos
        END DO
    END DO
    CALL bld_mtrstr(falt,name,nvbls,vrbl,units,rgn,info,dmngmtry,dataval,mtrstr)                                                    Finalmente construye el Metaraster
    CLOSE(ud)                                                                                                                       Cierra el archivo de datos
    CALL CPU_TIME(timeb)                                                                                                            Pregunta por la hora para estimar el tiempo de ejecuci�n del programa
    exetime=timeb-timea                                                                                                             Calcula el tiempo de ejecuci�n
    hrs=FLOOR(exetime/3600);mnts=FLOOR((exetime/3600-hrs)*60);scnds=exetime-mnts*60-hrs*3600                                        Convierte las unidades del tiempo de ejecuci�n
    WRITE(*,'(A101,2(I2,A2),F4.1,A6)')'*---|Importing Metaraster (*.metaRST, *.dataRST) files has been finished.'//&                Imprime el mensaje de finalizaci�n satisfactoria del programa
        & '            Execution time: ',hrs,'h ',mnts,'m ',scnds,'s|---*'
    !
END SUBROUTINE bld_mtrstrff_sp
!BL
SUBROUTINE bld_mtrstrff_dp(prjctrt,flnm,vrbl,units,X,Y,Z,dataval,mtrstr)
    !
    USE festellustype
    USE util,   ONLY : error
    USE iofrmtng, ONLY : strtonum
    USE domainprop, ONLY : bld_domainff
    USE raster, ONLY :  bld_mtrstr
    !
    IMPLICIT NONE
    !
    CHARACTER(*), INTENT(IN)                                        :: prjctrt,flnm
    CHARACTER(10), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(OUT)   :: units
    CHARACTER(20), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(OUT)   :: vrbl
    REAL(DP), DIMENSION(:,:), ALLOCATABLE, TARGET, INTENT(OUT)      :: dataval
    REAL(DP), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(OUT)        :: X,Y,Z
    TYPE(metaraster), INTENT(OUT)                                   :: mtrstr
    !
    CHARACTER(1)                                        :: enter
    CHARACTER(2)                                        :: rtdlmtng
    CHARACTER(20)                                       :: name,rgn,strnum
    CHARACTER(200)                                      :: rtm,rtd
    CHARACTER(500)                                      :: info
    INTEGER(I4B)                                        :: um,ud,i,j,nvbls,hrs,mnts
    REAL(DP)                                            :: timea,timeb,exetime,scnds,falt
    LOGICAL(LGT)                                        :: flexist
    TYPE(domain)                                        :: dmngmtry
    !
    PARAMETER(um=30,ud=31)
    !
    CALL CPU_TIME(timea)                                                                                                            Se pregunta por la hora para controlar el tiempo de ejecuci�n
    WRITE(*,'(A59)')'*---|Importing Metaraster (*.metaRST, *.dataRST) files|---*'                                                   Impresi�n del aviso de inicio de la importaci�n
    !Definici�n del caracter de delimitaci�n de rutas de acuerdo con el sistema operativo
    IF (INDEX(prjctrt,'/')==0) THEN
        rtdlmtng='\\'                                                                                                               Delimitador en MS-DOS
    ELSE
        rtdlmtng='//'                                                                                                               Delimitador en Unix-OS
    END IF
    !Verificaci�n de la existencia y pertinencia de los datos de entrada
    IF (LEN_TRIM(TRIM(ADJUSTL(prjctrt))//rtdlmtng//'Data'//rtdlmtng//'RST'//rtdlmtng//TRIM(ADJUSTL(flnm))//'.metaRST') > 200) THEN
        CALL error('in bld_mtrstrff_dp: the path of the Metaraster is very long')    
    ELSE
        rtm=TRIM(ADJUSTL(prjctrt))//rtdlmtng//'Data'//rtdlmtng//'RST'//rtdlmtng//TRIM(ADJUSTL(flnm))//'.metaRST'                    Archivo en el que  se almacenan los metadatos
        rtd=TRIM(ADJUSTL(prjctrt))//rtdlmtng//'Data'//rtdlmtng//'RST'//rtdlmtng//TRIM(ADJUSTL(flnm))//'.dataRST'                    Archivo en el que  se almacena el arreglo de datos
    END IF
    INQUIRE(FILE=TRIM(ADJUSTL(rtm)),EXIST=flexist)
    IF (.NOT.flexist) CALL error('in bld_mtrstrff_dp: Metaraster metadata file does not exist')                                     Verifica que el archivo con los metadatos exista y est� en la ruta correcta
    INQUIRE(FILE=TRIM(ADJUSTL(rtd)),EXIST=flexist)
    IF (.NOT.flexist) CALL error('in bld_mtrstrff_dp: Metaraster data file does not exist')                                         Verifica que el archivo con el arreglo de datos exista y est� en la ruta correcta
    !Lectura del archivo de metadatos
    OPEN(um,FILE=rtm,STATUS='OLD')                                                                                                  Abre el archivo de los metadatos
    READ(um,*);READ(um,'(A20)')name                                                                                                 Nombre del metaraster                                                                                    Entero que indica el tipo de n�mero para manejar los datos
    READ(um,*);READ(um,*);READ(um,'(A20)')strnum                                                                                  
    CALL strtonum(strnum,falt)                                                                                                      Valor de datos faltantes en precisi�n real doble
    READ(um,*);READ(um,*);READ(um,'(I9)')nvbls                                                                                      N�mero de variables
    ALLOCATE(vrbl(nvbls),units(nvbls))                                                                                              Aloja los arreglos de variables y unidades
    READ(um,*);READ(um,*)   
    DO i=1,nvbls
        READ(um,'(A20)')vrbl(i)                                                                                                     Nombres de la variables almacenadas en el Metaraster
    END DO
    READ(um,*);READ(um,*)
    DO i=1,nvbls
        READ(um,'(A10)')units(i)                                                                                                    Unidades de las variables almacenadas en el Metaraster
    END DO
    READ(um,*);READ(um,*);READ(um,'(A20)')rgn                                                                                       Regi�n en la que se define el Metaraster
    READ(um,*);READ(um,*);READ(um,'(A500)')info                                                                                     Informaci�n de inter�s sobre el Metaraster
    CLOSE(um)                                                                                                                       Cierra el archivo de los metadatos
    !Construye el dominio espacial
    CALL bld_domainff(prjctrt,flnm,X,Y,Z,dmngmtry)
    !Lectura del archivo que almacena el arreglo de datos
    ALLOCATE(dataval(dmngmtry%clmns*dmngmtry%rws*dmngmtry%lvls,nvbls))                                                              Aloja el arreglo de datos con precisi�n real doble
    OPEN(ud,FILE=rtd,STATUS='OLD')                                                                                                  Abre el archivo de datos
    DO j=1,nvbls
        DO i=1,dmngmtry%clmns*dmngmtry%rws*dmngmtry%lvls
            READ(ud,'(SP,ES18.11E2,A1)')dataval(i,j),enter                                                                          Y lo llena con la informaci�n del archivo de datos
        END DO
    END DO
    CALL bld_mtrstr(falt,name,nvbls,vrbl,units,rgn,info,dmngmtry,dataval,mtrstr)                                                    Finalmente construye el Metaraster
    CLOSE(ud)                                                                                                                       Cierra el archivo de datos
    CALL CPU_TIME(timeb)                                                                                                            Pregunta por la hora para estimar el tiempo de ejecuci�n del programa
    exetime=timeb-timea                                                                                                             Calcula el tiempo de ejecuci�n
    hrs=FLOOR(exetime/3600);mnts=FLOOR((exetime/3600-hrs)*60);scnds=exetime-mnts*60-hrs*3600                                        Convierte las unidades del tiempo de ejecuci�n
    WRITE(*,'(A101,2(I2,A2),F4.1,A6)')'*---|Importing Metaraster (*.metaRST, *.dataRST) files has been finished.'//&                Imprime el mensaje de finalizaci�n satisfactoria del programa
        & '            Execution time: ',hrs,'h ',mnts,'m ',scnds,'s|---*'
    !
END SUBROUTINE bld_mtrstrff_dp
!
!***********************************************************************************************************************************
!*********************** CREACI�N DE ACCESO A UNA VARIABLE TIPO METARASTER ALMACENADA ARCHIVO **************************************
!******************************** @authores �lvarez-Villa O.D. & Rend�n-�lvarez J.P. ***********************************************
!***********************************************************************************************************************************
!
!DESCRIPCI�N:
!Una subrutina que permite llenar los atributos de una variable tipo Metaraster que permiten el acceso a los archivos de texto
!(*.metaRST, *.domnRST y *.dataRST) almacenado en el disco duro
!
!DATOS DE ENTRADA:
!prjctrt        : Ruta del proyecto de trabajo                                                                            {caracter}
!flnm           : Nombre del proyecto de trabajo                                                                          {caracter}
!typdmn         : Entero que indica el tipo de n�mero para representar el dominio                                           {entero}
!typnum         : Entero que indica el tipo de n�mero para representar los datos                                            {entero}
!domflind       : Localizaci�n en memoria para localizar al archivo del dominio                                             {entero}
!datflind       : Localizaci�n en memoria para localizar al archivo de datos                                                {entero}
!
!DATOS DE SALIDA:
!mtrstr         : Variable tipo Metaraster a la cual se le asignan los par�metros                                       {metaraster}
!
SUBROUTINE bld_mtrstrif(prjctrt,flnm,typdmn,typnum,domflind,datflind,mtrstr)
    !
    USE festellustype
    USE util, ONLY: error
    USE domainprop, ONLY: bld_domainif
    !
    IMPLICIT NONE
    !
    CHARACTER(*), INTENT(IN)      :: prjctrt,flnm
    INTEGER(I4B), INTENT(IN)      :: typdmn,typnum,domflind,datflind
    TYPE(metaraster), INTENT(OUT) :: mtrstr
    !
    CHARACTER(2)   :: rtdlmtng
    CHARACTER(8)   :: header
    CHARACTER(200) :: rtd,rtm       
    INTEGER(I4B)   :: i,um,hrs,mnts
    REAL(DP)       :: timea,timeb,exetime,scnds
    LOGICAL(LGT)   :: flexist
    !
    PARAMETER(um=30)
    !
    CALL CPU_TIME(timea)                                                                                                            Se pregunta por la hora para controlar el tiempo de ejecuci�n
    WRITE(*,'(A75)')'*---|Starts creating access to Metaraster (*.metaRST, *.dataRST) files|---*'                                   Impresi�n del aviso de inicio de la importaci�n
    !Definici�n del caracter de delimitaci�n de rutas de acuerdo con el sistema operativo
    IF (INDEX(prjctrt,'/')==0) THEN
        rtdlmtng='\\'                                                                                                               Delimitador en MS-DOS
    ELSE
        rtdlmtng='//'                                                                                                               Delimitador en Unix-OS
    END IF
    !Verificaci�n de la existencia y pertinencia de los datos de entrada
    IF (LEN_TRIM(TRIM(ADJUSTL(prjctrt))//rtdlmtng//'Data'//rtdlmtng//'RST'//rtdlmtng//TRIM(ADJUSTL(flnm))//'.metaRST') > 200) THEN
        CALL error('in bld_mtrstrif: the path of the Metaraster is very long')    
    ELSE
        rtm=TRIM(ADJUSTL(prjctrt))//rtdlmtng//'Data'//rtdlmtng//'RST'//rtdlmtng//TRIM(ADJUSTL(flnm))//'.metaRST'                    Archivo en el que  se almacenan los metadatos
        rtd=TRIM(ADJUSTL(prjctrt))//rtdlmtng//'Data'//rtdlmtng//'RST'//rtdlmtng//TRIM(ADJUSTL(flnm))//'.dataRST'                    Archivo en el que  se almacena el arreglo de datos
    END IF
    INQUIRE(FILE=TRIM(ADJUSTL(rtm)),EXIST=flexist)
    IF (.NOT.flexist) CALL error('in bld_mtrstrif: Metaraster metadata file does not exist')                                        Verifica que el archivo con los metadatos exista y est� en la ruta correcta
    INQUIRE(FILE=TRIM(ADJUSTL(rtd)),EXIST=flexist)
    IF (.NOT.flexist) CALL error('in bld_mtrstrif: Metaraster data file does not exist')                                            Verifica que el archivo con el arreglo de datos exista y est� en la ruta correcta
    !Construcci�n del acceso al archivo del dominio
    CALL bld_domainif(prjctrt,flnm,domflind,mtrstr%dom)
    !Lectura del archivo de metadatos
    OPEN(um,FILE=rtm,STATUS='OLD',ACTION='READ')                                                                                    Abre el archivo de los metadatos
    READ(um,*);READ(um,'(A20)')mtrstr%name                                                                                          Nombre del metaraster                                                                                    Entero que indica el tipo de n�mero para manejar los datos
    SELECT CASE(typnum)                                                                                                             De acuerdo con el tipo de n�mero define
        CASE(0)
            READ(um,*);READ(um,*);READ(um,'(I9)')mtrstr%falti                                                                       Valor de datos faltantes en precisi�n entera
        CASE(1)
            READ(um,*);READ(um,*);READ(um,'(F16.8)')mtrstr%falts                                                                    Valor de datos faltantes en precisi�n real simple
        CASE(2)
            READ(um,*);READ(um,*);READ(um,'(F16.8)')mtrstr%faltd                                                                    Valor de datos faltantes en precisi�n real doble
        CASE DEFAULT
            CALL error('in bld_mtrstrif: non supported data number format')
    END SELECT
    READ(um,*);READ(um,*);READ(um,'(I9)')mtrstr%nvbles                                                                              N�mero de variables
    ALLOCATE(mtrstr%vrbl(mtrstr%nvbles),mtrstr%units(mtrstr%nvbles))                                                                Aloja los arreglos de variables y unidades
    READ(um,*);READ(um,*)   
    DO i=1,mtrstr%nvbles
        READ(um,'(A20)')mtrstr%vrbl(i)                                                                                              Nombres de la variables almacenadas en el Metaraster
    END DO
    READ(um,*);READ(um,*)
    DO i=1,mtrstr%nvbles
        READ(um,'(A10)')mtrstr%units(i)                                                                                             Unidades de las variables almacenadas en el Metaraster
    END DO
    READ(um,*);READ(um,*);READ(um,'(A20)')mtrstr%rgn                                                                                Regi�n en la que se define el Metaraster
    READ(um,*);READ(um,*);READ(um,'(A500)')mtrstr%info                                                                              Informaci�n de inter�s sobre el Metaraster
    CLOSE(um)                                                                                                                       Cierra el archivo de los metadatos
    !Llenado de los atributos para acceso a los archivos de metadatos y datos
    mtrstr%exist=.true.
    mtrstr%infile=.true.
    mtrstr%rtmet=rtm
    mtrstr%rtdat=rtd
    mtrstr%flind=datflind
    mtrstr%typdmn=typdmn
    mtrstr%typnum=typnum
    NULLIFY(mtrstr%datavali,mtrstr%datavalsp,mtrstr%datavaldp)
    CALL CPU_TIME(timeb)                                                                                                            Pregunta por la hora para estimar el tiempo de ejecuci�n del programa
    exetime=timeb-timea                                                                                                             Calcula el tiempo de ejecuci�n
    hrs=FLOOR(exetime/3600);mnts=FLOOR((exetime/3600-hrs)*60);scnds=exetime-mnts*60-hrs*3600                                        Convierte las unidades del tiempo de ejecuci�n
    WRITE(*,'(A111,2(I2,A2),F4.1,A6)')'*---|Creating acces to Metaraster (*.metaRST y *.dataRST) files has been'//&                 Imprime el mensaje de finalizaci�n satisfactoria del programa
        & '             finished. Execution time: ',hrs,'h ',mnts,'m ',scnds,'s|---*'
    !
END SUBROUTINE bld_mtrstrif
!
!***********************************************************************************************************************************
!************** CONSTRUCCI�N DE UNA VARIABLE TIPO METARASTER EN LA RAM CON BASE EN VARIABLES ALMACENADAS EN LA RAM *****************
!******************************** @authores �lvarez-Villa O.D. & Rend�n-�lvarez J.P. ***********************************************
!***********************************************************************************************************************************
!
!DESCRIPCI�N:
!Un conjunto de subrutinas que llenan los atributos del objeto tipo metaraster, a partir de informaci�n previamente almacenada en la
!memoria RAM
!
!DATOS DE ENTRADA:
!falt           : Valor de datos faltantes                                                                           {entero � real}
!name           : Nombre del campo                                                                                        {caracter}
!nvbls          : N�mero de variables almacenadas en el campo                                                               {entero}
!vrbl           : Arreglo con los nombres de la variables                                                           {caracter}(:)[O]
!units          : Arreglo con las unidades de la variables                                                          {caracter}(:)[O]
!rgn            : Regi�n geogr�fica en la que se ubica el campo                                                           {caracter}
!info           : Informaci�n de inter�s acerca del campo                                                                 {caracter}
!domn           : Variable tipo del dominio espacial                                                                        {domain}
!dataval        : Arreglo de valores para cuando el metaraster se almacena en RAM                              {entero � real}(:)[O]
!
!DATOS DE SALIDA:
!mtrstr         : Variable tipo metaraster

!
SUBROUTINE bld_mtrstr_i(falt,name,nvbls,vrbl,units,rgn,info,domn,dataval,mtrstr)    
    USE festellustype
    USE util, ONLY : error, assert_eq
    !
    IMPLICIT NONE
    !
    CHARACTER(20), INTENT(IN)                                       :: name,rgn
    CHARACTER(10), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(IN)    :: units
    CHARACTER(20), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(IN)    :: vrbl
    CHARACTER(500), INTENT(IN)                                      :: info
    INTEGER(I4B), INTENT(IN)                                        :: falt,nvbls
    INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE, TARGET, INTENT(IN)   :: dataval
    TYPE(domain), INTENT(IN)                                        :: domn
    TYPE(metaraster), INTENT(OUT)                                   :: mtrstr
    !
    INTEGER(I4B) :: n,nv,cnt
    !
    !Verificaci�n de consistencia en los datos de entrada
    IF ((SIZE(dataval,1)/=domn%clmns*domn%rws*domn%lvls).OR.(SIZE(dataval,2)/=nvbls))&
        & CALL error('in bld_mtrstr_i: incongruent size of data array')
    IF (SIZE(vrbl,1)/=nvbls) CALL error('in bld_mtrstr_i: incongruent size of variables names array')
    IF (SIZE(units,1)/=nvbls) CALL error('in bld_mtrstr_i: incongruent size of variables units array')
    !Construcci�n del Metaraster
    mtrstr%exist=.true.
    mtrstr%infile=.false.
    mtrstr%rtdat=''
    mtrstr%rtmet=''
    mtrstr%flind=-9999
    mtrstr%typdmn=1
    mtrstr%typnum=0
    mtrstr%nvbles=nvbls
    mtrstr%name=name
    mtrstr%vrbl=>vrbl
    mtrstr%units=>units
    mtrstr%rgn=rgn
    mtrstr%info=info
    mtrstr%falts=falt
    mtrstr%dom%clmns=domn%clmns
    mtrstr%dom%rws=domn%rws
    mtrstr%dom%lvls=domn%lvls
    mtrstr%dom%exist=.true.
    IF (.NOT.domn%infile) THEN 
        mtrstr%dom%infile=.false.
        mtrstr%dom%rtdom=''
        mtrstr%dom%flind=-9999
        mtrstr%dom%Xs=>domn%Xs
        mtrstr%dom%Ys=>domn%Ys
        mtrstr%dom%Zs=>domn%Zs
    ELSE
        mtrstr%dom%infile=.true.
        mtrstr%dom%rtdom=domn%rtdom
        mtrstr%dom%flind=domn%flind
        NULLIFY(mtrstr%dom%Xs,mtrstr%dom%Ys,mtrstr%dom%Zs,mtrstr%dom%Xd,mtrstr%dom%Yd,mtrstr%dom%Zd)
    END IF
    mtrstr%datavali=>dataval
    NULLIFY(mtrstr%datavalsp,mtrstr%datavaldp)
    !
END SUBROUTINE bld_mtrstr_i
!BL
SUBROUTINE bld_mtrstr_sp(falt,name,nvbls,vrbl,units,rgn,info,domn,dataval,mtrstr)    
    USE festellustype
    USE util, ONLY : error, assert_eq
    !
    IMPLICIT NONE
    !
    CHARACTER(20), INTENT(IN)                                       :: name,rgn
    CHARACTER(10), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(IN)    :: units
    CHARACTER(20), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(IN)    :: vrbl
    CHARACTER(500), INTENT(IN)                                      :: info
    INTEGER(I4B), INTENT(IN)                                        :: nvbls
    REAL(SP), INTENT(IN)                                            :: falt
    REAL(SP), DIMENSION(:,:), ALLOCATABLE, TARGET, INTENT(IN)       :: dataval
    TYPE(domain), INTENT(IN)                                        :: domn
    TYPE(metaraster), INTENT(OUT)                                   :: mtrstr
    !
    INTEGER(I4B) :: n,nv,cnt
    !
    !Verificaci�n de consistencia en los datos de entrada
    IF ((SIZE(dataval,1)/=domn%clmns*domn%rws*domn%lvls).OR.(SIZE(dataval,2)/=nvbls))&
        & CALL error('in bld_mtrstr_sp: incongruent size of data array')
    IF (SIZE(vrbl,1)/=nvbls) CALL error('in bld_mtrstr_sp: incongruent size of variables names array')
    IF (SIZE(units,1)/=nvbls) CALL error('in bld_mtrstr_sp: incongruent size of variables units array')
    !Construcci�n del Metaraster
    mtrstr%exist=.true.
    mtrstr%infile=.false.
    mtrstr%rtdat=''
    mtrstr%rtmet=''
    mtrstr%flind=-9999
    mtrstr%typdmn=1
    mtrstr%typnum=1
    mtrstr%nvbles=nvbls
    mtrstr%name=name
    mtrstr%vrbl=>vrbl
    mtrstr%units=>units
    mtrstr%rgn=rgn
    mtrstr%info=info
    mtrstr%falts=falt
    mtrstr%dom%clmns=domn%clmns
    mtrstr%dom%rws=domn%rws
    mtrstr%dom%lvls=domn%lvls
    mtrstr%dom%exist=.true.
    IF (.NOT.domn%infile) THEN 
        mtrstr%dom%infile=.false.
        mtrstr%dom%rtdom=''
        mtrstr%dom%flind=-9999
        mtrstr%dom%Xs=>domn%Xs
        mtrstr%dom%Ys=>domn%Ys
        mtrstr%dom%Zs=>domn%Zs
    ELSE
        mtrstr%dom%infile=.true.
        mtrstr%dom%rtdom=domn%rtdom
        mtrstr%dom%flind=domn%flind
        NULLIFY(mtrstr%dom%Xs,mtrstr%dom%Ys,mtrstr%dom%Zs,mtrstr%dom%Xd,mtrstr%dom%Yd,mtrstr%dom%Zd)
    END IF
    mtrstr%datavalsp=>dataval
    NULLIFY(mtrstr%datavali,mtrstr%datavaldp)
    !
END SUBROUTINE bld_mtrstr_sp
!BL
SUBROUTINE bld_mtrstr_dp(falt,name,nvbls,vrbl,units,rgn,info,domn,dataval,mtrstr)    
    USE festellustype
    USE util, ONLY : error, assert_eq
    !
    IMPLICIT NONE
    !
    CHARACTER(20), INTENT(IN)                                       :: name,rgn
    CHARACTER(10), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(IN)    :: units
    CHARACTER(20), DIMENSION(:), ALLOCATABLE, TARGET, INTENT(IN)    :: vrbl
    CHARACTER(500), INTENT(IN)                                      :: info
    INTEGER(I4B), INTENT(IN)                                        :: nvbls
    REAL(DP), INTENT(IN)                                            :: falt
    REAL(DP), DIMENSION(:,:), ALLOCATABLE, TARGET, INTENT(IN)       :: dataval
    TYPE(domain), INTENT(IN)                                        :: domn
    TYPE(metaraster), INTENT(OUT)                                   :: mtrstr
    !
    INTEGER(I4B) :: n,nv,cnt
    !
    !Verificaci�n de consistencia en los datos de entrada
    IF ((SIZE(dataval,1)/=domn%clmns*domn%rws*domn%lvls).OR.(SIZE(dataval,2)/=nvbls))&
        & CALL error('in bld_mtrstr_dp: incongruent size of data array')
    IF (SIZE(vrbl,1)/=nvbls) CALL error('in bld_mtrstr_dp: incongruent size of variables names array')
    IF (SIZE(units,1)/=nvbls) CALL error('in bld_mtrstr_dp: incongruent size of variables units array')
    !Construcci�n del Metaraster
    mtrstr%exist=.true.
    mtrstr%infile=.false.
    mtrstr%rtdat=''
    mtrstr%rtmet=''
    mtrstr%flind=-9999
    mtrstr%typdmn=2
    mtrstr%typnum=2
    mtrstr%nvbles=nvbls
    mtrstr%name=name
    mtrstr%vrbl=>vrbl
    mtrstr%units=>units
    mtrstr%rgn=rgn
    mtrstr%info=info
    mtrstr%faltd=falt
    mtrstr%dom%clmns=domn%clmns
    mtrstr%dom%rws=domn%rws
    mtrstr%dom%lvls=domn%lvls
    mtrstr%dom%exist=.true.
    IF (.NOT.domn%infile) THEN 
        mtrstr%dom%infile=.false.
        mtrstr%dom%rtdom=''
        mtrstr%dom%flind=-9999
        mtrstr%dom%Xd=>domn%Xd
        mtrstr%dom%Yd=>domn%Yd
        mtrstr%dom%Zd=>domn%Zd
    ELSE
        mtrstr%dom%infile=.true.
        mtrstr%dom%rtdom=domn%rtdom
        mtrstr%dom%flind=domn%flind
        NULLIFY(mtrstr%dom%Xs,mtrstr%dom%Ys,mtrstr%dom%Zs,mtrstr%dom%Xd,mtrstr%dom%Yd,mtrstr%dom%Zd)
    END IF
    mtrstr%datavaldp=>dataval
    NULLIFY(mtrstr%datavali,mtrstr%datavalsp)
    !
END SUBROUTINE bld_mtrstr_dp
!
!***********************************************************************************************************************************
!*********************** EXPORTACI�N DE UNA VARIABLE TIPO METARASTER A ARCHIVOS EN EL DISCO DURO ***********************************
!******************************* @authores: Alvarez-Villa O.D. & Rend�n-�lvarez J.P. ***********************************************
!***********************************************************************************************************************************
!
!DESCRIPCI�N:
!Una subrutina que permite exportar la informaci�n de una variable tipo metaraster que se encuentra en la RAM a archivos de texto
!(*.metaRST,*.domnRST, *.dataRST) en el disco duro
!
!DATOS DE ENTRADA:
!prjctrt    : Ruta del proyecto de trabajo                                                                                {caracter}
!fldr       : Carpeta del proyecto en la que almacenar� el Metaraster                                                     {caracter}
!flnm       : Nombre con el que se almacenar� el archivo en el disco duro                                                 {caracter}
!mtrstr     : Variable tipo metaraster a exportar                                                                       {metaraster}
!
SUBROUTINE svmetaraster(prjctrt,fldr,flnm,mtrstr)
    USE festellustype
    USE util, ONLY: error
    USE tpvrblchckng, ONLY: chckmtrstr
    USE iofrmtng, ONLY: bldwrtfrmt,bldflnm
    USE domainprop, ONLY: svdomain
    !
    IMPLICIT NONE
    !
    CHARACTER(*), INTENT(IN)         :: prjctrt,fldr
    CHARACTER(*), INTENT(INOUT)      :: flnm
    TYPE(metaraster), INTENT(INOUT)  :: mtrstr
    !
    CHARACTER(2)    :: rtdlmtng
    CHARACTER(8)    :: hdfrmt
    CHARACTER(200)  :: rtm,rtd
    INTEGER(I4B)    :: i,j,k,um,ud,hrs,mnts
    REAL(DP)        :: timea,timeb,exetime,scnds
    !
    PARAMETER(um=88,ud=89)
    !
    CALL CPU_TIME(timea)                                                                                                            Se pregunta por la hora para controlar el tiempo de ejecuci�n
    WRITE(*,'(A62)')'*---|Exporting Metaraster to (*.metaRST, *.dataRST) files|---*'                                                Impresi�n del aviso de inicio de la importaci�n
    !Definici�n del caracter de delimitaci�n de rutas de acuerdo con el sistema operativo
    IF (INDEX(prjctrt,'/')==0) THEN
        rtdlmtng='\\'                                                                                                               Delimitador en MS-DOS
    ELSE
        rtdlmtng='//'                                                                                                               Delimitador en Unix-OS
    END IF
    !Verificaci�n de la consistencia en los datos de entrada
    IF ((fldr/='Data').AND.(fldr/='Results')) CALL error('in svmetaraster: non allowed storage folder')
    IF (LEN_TRIM(TRIM(ADJUSTL(prjctrt))//rtdlmtng//fldr//rtdlmtng//'RST'//rtdlmtng//TRIM(ADJUSTL(flnm))//'.metaRST')>196) THEN 
        CALL error('in svmetaraster: the path of te Metaraster is very long')
    ELSE
        CALL bldflnm(TRIM(ADJUSTL(prjctrt))//rtdlmtng//fldr//rtdlmtng//'RST',flnm,'metaRST')
        CALL bldflnm(TRIM(ADJUSTL(prjctrt))//rtdlmtng//fldr//rtdlmtng//'RST',flnm,'dataRST')
        CALL bldflnm(TRIM(ADJUSTL(prjctrt))//rtdlmtng//fldr//rtdlmtng//'RST',flnm,'domnRST')
        rtm=TRIM(ADJUSTL(prjctrt))//rtdlmtng//fldr//rtdlmtng//'RST'//rtdlmtng//TRIM(ADJUSTL(flnm))//'.metaRST'
        rtd=TRIM(ADJUSTL(prjctrt))//rtdlmtng//fldr//rtdlmtng//'RST'//rtdlmtng//TRIM(ADJUSTL(flnm))//'.dataRST'
    END IF
    CALL chckmtrstr(mtrstr,mtrstr%typnum,'svmetaraster')
    !Escritura del archivo de metadatos (*.metaRST)
    OPEN(um,FILE=rtm,STATUS='NEW')
    WRITE(um,'(A17)')'[METARASTER NAME]'
    CALL bldwrtfrmt(mtrstr%name,hdfrmt);WRITE(um,'('//TRIM(ADJUSTL(hdfrmt))//')')TRIM(ADJUSTL(mtrstr%name));WRITE(um,*)
    SELECT CASE(mtrstr%typnum)
        CASE(0)
            WRITE(um,'(A14)')'[NODATA VALUE]'
            CALL bldwrtfrmt(mtrstr%falti,hdfrmt);WRITE(um,'('//TRIM(ADJUSTL(hdfrmt))//')')mtrstr%falti;WRITE(um,*)
        CASE(1)
            WRITE(um,'(A14)')'[NODATA VALUE]'
            CALL bldwrtfrmt(mtrstr%falts,hdfrmt);WRITE(um,'('//TRIM(ADJUSTL(hdfrmt))//')')mtrstr%falts;WRITE(um,*)
        CASE(2)
            WRITE(um,'(A14)')'[NODATA VALUE]'
            CALL bldwrtfrmt(mtrstr%faltd,hdfrmt);WRITE(um,'('//TRIM(ADJUSTL(hdfrmt))//')')mtrstr%faltd;WRITE(um,*)
        CASE DEFAULT
            CALL error('in svmetaraster: non supported data value')
    END SELECT
    WRITE(um,'(A21)')'[NUMBER OF VARIABLES]'
    CALL bldwrtfrmt(mtrstr%nvbles,hdfrmt);WRITE(um,'('//TRIM(ADJUSTL(hdfrmt))//')')mtrstr%nvbles;WRITE(um,*)
    WRITE(um,'(A17)')'[VARIABLES NAMES]'
    DO i=1,mtrstr%nvbles
        CALL bldwrtfrmt(mtrstr%vrbl(i),hdfrmt);WRITE(um,'('//TRIM(ADJUSTL(hdfrmt))//')')TRIM(ADJUSTL(mtrstr%vrbl(i)))
    END DO
    WRITE(um,*);WRITE(um,'(A17)')'[VARIABLES UNITS]'
    DO i=1,mtrstr%nvbles
        CALL bldwrtfrmt(mtrstr%units(i),hdfrmt);WRITE(um,'('//TRIM(ADJUSTL(hdfrmt))//')')TRIM(ADJUSTL(mtrstr%units(i)))
    END DO
    WRITE(um,*);WRITE(um,'(A8)')'[REGION]'
    CALL bldwrtfrmt(mtrstr%rgn,hdfrmt);WRITE(um,'('//TRIM(ADJUSTL(hdfrmt))//')')TRIM(ADJUSTL(mtrstr%rgn));WRITE(um,*)
    WRITE(um,'(A13)')'[INFORMATION]'
    CALL bldwrtfrmt(mtrstr%info,hdfrmt);WRITE(um,'('//TRIM(ADJUSTL(hdfrmt))//')')TRIM(ADJUSTL(mtrstr%info))
    CLOSE(um)
    !Escritura del archivo de domino (*.domnRST)
    CALL svdomain(prjctrt,fldr,flnm,mtrstr%dom,mtrstr%typdmn)
    !Escritura del archivo de datos (*.dataRST)
    SELECT CASE(mtrstr%typnum)
        CASE(0)
            OPEN(ud,FILE=rtd,STATUS='NEW',ACCESS='DIRECT', RECL=12, FORM='FORMATTED')
            DO j=1,mtrstr%nvbles
                DO i=1,mtrstr%dom%clmns*mtrstr%dom%rws*mtrstr%dom%lvls
                    k=(j-1)*mtrstr%dom%clmns*mtrstr%dom%rws*mtrstr%dom%lvls+i
                    WRITE(ud,'(SP,I11,A1)',REC=k)mtrstr%datavali(i,j),13
                END DO
            END DO
        CASE(1)
            OPEN(ud,FILE=rtd,STATUS='NEW',ACCESS='DIRECT', RECL=19, FORM='FORMATTED')
            DO j=1,mtrstr%nvbles
                DO i=1,mtrstr%dom%clmns*mtrstr%dom%rws*mtrstr%dom%lvls
                    k=(j-1)*mtrstr%dom%clmns*mtrstr%dom%rws*mtrstr%dom%lvls+i
                    WRITE(ud,'(SP,ES18.11E2,A1)',REC=k)mtrstr%datavalsp(i,j),13
                END DO
            END DO
        CASE(2)
            OPEN(ud,FILE=rtd,STATUS='NEW',ACCESS='DIRECT', RECL=19, FORM='FORMATTED')
            DO j=1,mtrstr%nvbles
                DO i=1,mtrstr%dom%clmns*mtrstr%dom%rws*mtrstr%dom%lvls
                    k=(j-1)*mtrstr%dom%clmns*mtrstr%dom%rws*mtrstr%dom%lvls+i
                    WRITE(ud,'(SP,ES18.11E2,A1)',REC=k)mtrstr%datavaldp(i,j),13
                END DO
            END DO
        CASE DEFAULT
            CALL error('in svmetaraster: non supported data value')
    END SELECT
    CLOSE(ud)
    CALL CPU_TIME(timeb)                                                                                                            Pregunta por la hora para estimar el tiempo de ejecuci�n del programa
    exetime=timeb-timea                                                                                                             Calcula el tiempo de ejecuci�n
    hrs=FLOOR(exetime/3600);mnts=FLOOR((exetime/3600-hrs)*60);scnds=exetime-mnts*60-hrs*3600                                        Convierte las unidades del tiempo de ejecuci�n
    WRITE(*,'(A116,2(I2,A2),F4.1,A6)')'*---|Exporting Metaraster to (*.metaRST, *.domnRST y *.dataRST) files has'//&                Imprime el mensaje de finalizaci�n satisfactoria del programa
        & '            been finished. Execution time: ',hrs,'h ',mnts,'m ',scnds,'s|---*'
    !
END SUBROUTINE svmetaraster
!
!***********************************************************************************************************************************
!*********************************** DESTRUCCION DE LA VARIABLE TIPO METARASTER ****************************************************
!******************************* @authores: �lvarez-Villa O.D. & Rend�n-�lvarez J.P. ***********************************************
!***********************************************************************************************************************************
!
!DESCRIPCI�N:
!Una rutina que desarma una variable tipo Metaraster. Cuando el Metaraster se haya construido en la memoria RAM lo que hace �ste
!algoritmo es llenar con los valores de faltantes los atributos que corresponden a metadatos y romper los punteros de los arreglos
!de datos y de coordenadas. Es importante aclarar que los arreglos no se desalojan, y en consecuencia la liberaci�n de memoria RAM 
!no es significativa; para liberar la RAM que ocupa el Metaraster deben desalojarse los objetivos a los que apuntan los arreglos de 
!datos y de coordenadas. Cuando el Metaraster se haya construido en el disco duro, pueden o no se rompen los v�nculos que 
!proporcionan acceso a los archivos con los arreglos de datos y de coordenadas
!
!DATOS DE ENTRADA:
!mtrstr         : Variable tipo Metaraster a destruir                                                                   {metaraster}
!switch         : Cuando es verdadera se destruye el Metaraster, bien sea que se haya construido en la RAM                  {l�gica}
!                 o en el disco duro. Cuando es falsa s�lo se destruye en la RAM
!
SUBROUTINE dstry_metaraster(mtrstr,switch)   
    !
    USE festellustype
    USE domainprop, ONLY : dstry_domain
    !
    IMPLICIT NONE
    !
    LOGICAL(LGT), INTENT(IN)        :: switch
    TYPE(metaraster), INTENT(INOUT) :: mtrstr
    !
    IF (switch) THEN
        mtrstr%exist=.false.
        mtrstr%infile=.false.
        mtrstr%rtdat=''
        mtrstr%rtmet=''
        mtrstr%flind=-9999
        mtrstr%typdmn=-9999   
        mtrstr%typnum=-9999       
        mtrstr%name=''
        mtrstr%nvbles=-9999
        NULLIFY(mtrstr%vrbl,mtrstr%units)
        mtrstr%rgn=''
        mtrstr%info=''
        mtrstr%falti=-1E9
        mtrstr%falts=-1.E9_sp
        mtrstr%faltd=-1.E9_dp
    END IF
    CALL dstry_domain(mtrstr%dom,switch)
    IF(ASSOCIATED(mtrstr%datavali))THEN
        NULLIFY(mtrstr%datavali)
    END IF
    IF(ASSOCIATED(mtrstr%datavalsp))THEN
        NULLIFY(mtrstr%datavalsp)
    END IF
    IF(ASSOCIATED(mtrstr%datavaldp))THEN
        NULLIFY(mtrstr%datavaldp)
    END IF
    !
END SUBROUTINE dstry_metaraster
!
!***********************************************************************************************************************************
!*************************** LECTURA DE UN VALOR EN LOS ARREGLOS DE DATOS DE UN METARASTER *****************************************
!********************************* @authores: �lvarez-Villa O.D. & Rend�n-�lvarez J.P.  ********************************************
!***********************************************************************************************************************************
!
!DESCRIPCI�N:
!Un conjunto de subrutinas que retornan un valor de un Metaraster almacenado en la RAM o en un archivo (*.dataRST) en el
!disco duro, con base en sus �ndices correspondientes de columna, fila, nivel y variable, o en el �ndice de almacenamiento como
!un vector columna
!
!DATOS DE ENTRADA:
!mtrstr             : Variable tipo Metaraster construida en archivo                                                    {metaraster}
!i                  : Columna del v�xel donde se encuentra el dato                                                [opcional]{entero}
!j                  : Fila del v�xel donde se encuentra el dato                                                   [opcional]{entero}
!k                  : Nivel del v�xel donde se encuentra el dato                                                  [opcional]{entero}
!v                  : �ndice de la variable a leer en el arreglo de datos                                         [opcional]{entero}
!l                  : �ndice integrado de almacenamiento como un vector columna                                   [opcional]{entero}        
!switch             : Variable l�gica que sirve para indicar si se desea revisar la adecuada construcci�n del Metaraster    {l�gico}
!
!DATOS DE SALIDA:
!vl                 : Valor le�do                                                                                    {entero � real}
!
SUBROUTINE gt_mtrstrvl_i(mtrstr,switch,vl,i,j,k,v,l)
    USE festellustype
    USE tpvrblchckng, ONLY : chckmtrstr
    USE util, ONLY : error
    !
    IMPLICIT NONE        
    !
    INTEGER(I4B), INTENT(OUT)              :: vl
    INTEGER(I4B), OPTIONAL, INTENT(IN)     :: i,j,k,v,l
    LOGICAL(LGT)                           :: switch
    TYPE(metaraster), INTENT(IN)           :: mtrstr
    !
    CHARACTER(1) :: enter
    INTEGER(I4B) :: pstn=0    
    !
    !Verificaci�n de la consistencia y pertinencia de los datos de entrada
    IF (switch) CALL chckmtrstr(mtrstr,0,'gt_mtrstrvl_i')                                                                           El Metaraster
    !Obtenci�n del valor
    IF (PRESENT(i).AND.PRESENT(j).AND.PRESENT(k).AND.PRESENT(v).AND.(.NOT.PRESENT(l)).AND.(mtrstr%infile.EQV.(.false.))) THEN       Cuando el Metaraster est� construido en la RAM y se especifica la columna, la fila y el nivel del v�xel y la variable
        pstn=(k-1)*mtrstr%dom%clmns*mtrstr%dom%rws+(j-1)*mtrstr%dom%clmns+i                                                         Calcula la posici�n del dato
        IF ((pstn<1).OR.(pstn>mtrstr%dom%clmns*mtrstr%dom%rws*mtrstr%dom%lvls).OR.(v<1).OR.(v>mtrstr%nvbles)) &                     Sale con error si los �ndices est�n fuera de rango
            & CALL error('gt_mtrstrvl_i: voxel index out of range')
        vl=mtrstr%datavali(pstn,v)                                                                                                  Valor a retornar
    ELSE IF (PRESENT(v).AND.PRESENT(l).AND.(.NOT.PRESENT(i)).AND.(.NOT.PRESENT(j))&                                                 Cuando el Metaraster est� construido en la RAM y se especifica el �ndice del v�xel en el vector columna y la variable 
            & .AND.(.NOT.PRESENT(k)).AND.(mtrstr%infile.EQV.(.false.))) THEN
        IF ((l<1).OR.(l>mtrstr%dom%clmns*mtrstr%dom%rws*mtrstr%dom%lvls).OR.(v<1).OR.(v>mtrstr%nvbles)) &                           Sale con error si los �ndices est�n fuera de rango
            & CALL error('gt_mtrstrvl_i: voxel index out of range')
        vl=mtrstr%datavali(l,v)                                                                                                     Valor a retornar
    ELSE IF (PRESENT(i).AND.PRESENT(j).AND.PRESENT(k).AND.PRESENT(v).AND.(.NOT.PRESENT(l)).AND.(mtrstr%infile.EQV.(.true.))) THEN   Cuando el Metaraster est� construido en el disco duro y se especifica la columna, la fila y el nivel del v�xel y la variable
        pstn=(v-1)*mtrstr%dom%clmns*mtrstr%dom%rws*mtrstr%dom%lvls+(k-1)*mtrstr%dom%clmns*mtrstr%dom%rws+(j-1)*mtrstr%dom%clmns+i   Calcula la posici�n del dato
        OPEN(mtrstr%flind,FILE=mtrstr%rtdat,ACTION='READ',STATUS='OLD',ACCESS='DIRECT',RECL=12,FORM='FORMATTED')                    Abre el archivo de datos
        READ(mtrstr%flind,'(SP,I11,A1)',REC=pstn)vl,enter                                                                           Leee el valor a retornar
        CLOSE(mtrstr%flind)                                                                                                         Cierra el archivo de datos
    ELSE IF (PRESENT(l).AND.(.NOT.PRESENT(v)).AND.(.NOT.PRESENT(i)).AND.(.NOT.PRESENT(j))&                                          Cuando el Metaraster est� construido en el disco duro y se especifica el �ndice del v�xel en el vector columna
            & .AND.(.NOT.PRESENT(k)).AND.(mtrstr%infile.EQV.(.true.))) THEN                                                 
        OPEN(mtrstr%flind,FILE=mtrstr%rtdat,ACTION='READ',STATUS='OLD',ACCESS='DIRECT',RECL=12,FORM='FORMATTED')                    Abre el archivo de datos                
        READ(mtrstr%flind,'(SP,I11,A1)',REC=l)vl,enter                                                                              Lee el valor a retornar
        CLOSE(mtrstr%flind)                                                                                                         Cierra el archivo de datos
    ELSE
        CALL error('gt_mtrstrvl_i: non supported indexes combination')                                                              Sale con error si no se soporta la combinaci�n de �ndices
    END IF   
    !
END SUBROUTINE gt_mtrstrvl_i
!BL
SUBROUTINE gt_mtrstrvl_sp(mtrstr,switch,vl,i,j,k,v,l)
    USE festellustype
    USE tpvrblchckng, ONLY : chckmtrstr
    USE util, ONLY : error
    !
    IMPLICIT NONE        
    !
    REAL(SP), INTENT(OUT)                  :: vl
    INTEGER(I4B), OPTIONAL, INTENT(IN)     :: i,j,k,v,l
    LOGICAL(LGT)                           :: switch
    TYPE(metaraster), INTENT(IN)           :: mtrstr
    !
    CHARACTER(1) :: enter
    INTEGER(I4B) :: pstn=0    
    !
    !Verificaci�n de la consistencia y pertinencia de los datos de entrada
    IF (switch) CALL chckmtrstr(mtrstr,1,'gt_mtrstrvl_sp')                                                                          El Metaraster
    !Obtenci�n del valor
    IF (PRESENT(i).AND.PRESENT(j).AND.PRESENT(k).AND.PRESENT(v).AND.(.NOT.PRESENT(l)).AND.(mtrstr%infile.EQV.(.false.))) THEN       Cuando el Metaraster est� construido en la RAM y se especifica la columna, la fila y el nivel del v�xel y la variable
        pstn=(k-1)*mtrstr%dom%clmns*mtrstr%dom%rws+(j-1)*mtrstr%dom%clmns+i                                                         Calcula la posici�n del dato
        IF ((pstn<1).OR.(pstn>mtrstr%dom%clmns*mtrstr%dom%rws*mtrstr%dom%lvls).OR.(v<1).OR.(v>mtrstr%nvbles)) &                     Sale con error si los �ndices est�n fuera de rango
            & CALL error('gt_mtrstrvl_sp: variable index out of range')
        vl=mtrstr%datavalsp(pstn,v)                                                                                                 Valor a retornar
    ELSE IF (PRESENT(v).AND.PRESENT(l).AND.(.NOT.PRESENT(i)).AND.(.NOT.PRESENT(j))&                                                 Cuando el Metaraster est� construido en la RAM y se especifica el �ndice del v�xel en el vector columna y la variable 
            & .AND.(.NOT.PRESENT(k)).AND.(mtrstr%infile.EQV.(.false.))) THEN
        IF ((l<1).OR.(l>mtrstr%dom%clmns*mtrstr%dom%rws*mtrstr%dom%lvls).OR.(v<1).OR.(v>mtrstr%nvbles)) &                           Sale con error si los �ndices est�n fuera de rango
            & CALL error('gt_mtrstrvl_sp: variable index out of range')
        vl=mtrstr%datavalsp(l,v)                                                                                                    Valor a retornar
    ELSE IF (PRESENT(i).AND.PRESENT(j).AND.PRESENT(k).AND.PRESENT(v).AND.(.NOT.PRESENT(l)).AND.(mtrstr%infile.EQV.(.true.))) THEN   Cuando el Metaraster est� construido en el disco duro y se especifica la columna, la fila y el nivel del v�xel y la variable
        pstn=(v-1)*mtrstr%dom%clmns*mtrstr%dom%rws*mtrstr%dom%lvls+(k-1)*mtrstr%dom%clmns*mtrstr%dom%rws+(j-1)*mtrstr%dom%clmns+i   Calcula la posici�n del dato
        OPEN(mtrstr%flind,FILE=mtrstr%rtdat,ACTION='READ',STATUS='OLD',ACCESS='DIRECT',RECL=19,FORM='FORMATTED')                    Abre el archivo de datos
        READ(mtrstr%flind,'(SP,ES18.11E2,A1)',REC=pstn)vl,enter                                                                     Leee el valor a retornar
        CLOSE(mtrstr%flind)                                                                                                         Cierra el archivo de datos
    ELSE IF (PRESENT(l).AND.(.NOT.PRESENT(v)).AND.(.NOT.PRESENT(i)).AND.(.NOT.PRESENT(j))&                                          Cuando el Metaraster est� construido en el disco duro y se especifica el �ndice del v�xel en el vector columna
            & .AND.(.NOT.PRESENT(k)).AND.(mtrstr%infile.EQV.(.true.))) THEN                                                 
        OPEN(mtrstr%flind,FILE=mtrstr%rtdat,ACTION='READ',STATUS='OLD',ACCESS='DIRECT',RECL=19,FORM='FORMATTED')                    Abre el archivo de datos                
        READ(mtrstr%flind,'(SP,ES18.11E2,A1)',REC=l)vl,enter                                                                        Lee el valor a retornar
        CLOSE(mtrstr%flind)                                                                                                         Cierra el archivo de datos
    ELSE
        CALL error('gt_mtrstrvl_sp: non supported indexes combination')                                                             Sale con error si no se soporta la combinaci�n de �ndices
    END IF   
    !
END SUBROUTINE gt_mtrstrvl_sp
!BL
SUBROUTINE gt_mtrstrvl_dp(mtrstr,switch,vl,i,j,k,v,l)
    USE festellustype
    USE tpvrblchckng, ONLY : chckmtrstr
    USE util, ONLY : error
    !
    IMPLICIT NONE        
    !
    REAL(DP), INTENT(OUT)                  :: vl
    INTEGER(I4B), OPTIONAL, INTENT(IN)     :: i,j,k,v,l
    LOGICAL(LGT)                           :: switch
    TYPE(metaraster), INTENT(IN)           :: mtrstr
    !
    CHARACTER(1) :: enter
    INTEGER(I4B) :: pstn=0    
    !
    !Verificaci�n de la consistencia y pertinencia de los datos de entrada
    IF (switch) CALL chckmtrstr(mtrstr,2,'gt_mtrstrvl_sp')                                                                          El Metaraster
    !Obtenci�n del valor
    IF (PRESENT(i).AND.PRESENT(j).AND.PRESENT(k).AND.PRESENT(v).AND.(.NOT.PRESENT(l)).AND.(mtrstr%infile.EQV.(.false.))) THEN       Cuando el Metaraster est� construido en la RAM y se especifica la columna, la fila y el nivel del v�xel y la variable
        pstn=(k-1)*mtrstr%dom%clmns*mtrstr%dom%rws+(j-1)*mtrstr%dom%clmns+i                                                         Calcula la posici�n del dato
        IF ((pstn<1).OR.(pstn>mtrstr%dom%clmns*mtrstr%dom%rws*mtrstr%dom%lvls).OR.(v<1).OR.(v>mtrstr%nvbles)) &                     Sale con error si los �ndices est�n fuera de rango
            & CALL error('gt_mtrstrvl_sp: variable index out of range')
        vl=mtrstr%datavaldp(pstn,v)                                                                                                 Valor a retornar
    ELSE IF (PRESENT(v).AND.PRESENT(l).AND.(.NOT.PRESENT(i)).AND.(.NOT.PRESENT(j))&                                                 Cuando el Metaraster est� construido en la RAM y se especifica el �ndice del v�xel en el vector columna y la variable 
            & .AND.(.NOT.PRESENT(k)).AND.(mtrstr%infile.EQV.(.false.))) THEN
        IF ((l<1).OR.(l>mtrstr%dom%clmns*mtrstr%dom%rws*mtrstr%dom%lvls).OR.(v<1).OR.(v>mtrstr%nvbles)) &                           Sale con error si los �ndices est�n fuera de rango
            & CALL error('gt_mtrstrvl_sp: variable index out of range')
        vl=mtrstr%datavaldp(l,v)                                                                                                    Valor a retornar
    ELSE IF (PRESENT(i).AND.PRESENT(j).AND.PRESENT(k).AND.PRESENT(v).AND.(.NOT.PRESENT(l)).AND.(mtrstr%infile.EQV.(.true.))) THEN   Cuando el Metaraster est� construido en el disco duro y se especifica la columna, la fila y el nivel del v�xel y la variable
        pstn=(v-1)*mtrstr%dom%clmns*mtrstr%dom%rws*mtrstr%dom%lvls+(k-1)*mtrstr%dom%clmns*mtrstr%dom%rws+(j-1)*mtrstr%dom%clmns+i   Calcula la posici�n del dato
        OPEN(mtrstr%flind,FILE=mtrstr%rtdat,ACTION='READ',STATUS='OLD',ACCESS='DIRECT',RECL=19,FORM='FORMATTED')                    Abre el archivo de datos
        READ(mtrstr%flind,'(SP,ES18.11E2,A1)',REC=pstn)vl,enter                                                                     Leee el valor a retornar
        CLOSE(mtrstr%flind)                                                                                                         Cierra el archivo de datos
    ELSE IF (PRESENT(l).AND.(.NOT.PRESENT(v)).AND.(.NOT.PRESENT(i)).AND.(.NOT.PRESENT(j))&                                          Cuando el Metaraster est� construido en el disco duro y se especifica el �ndice del v�xel en el vector columna
            & .AND.(.NOT.PRESENT(k)).AND.(mtrstr%infile.EQV.(.true.))) THEN                                                 
        OPEN(mtrstr%flind,FILE=mtrstr%rtdat,ACTION='READ',STATUS='OLD',ACCESS='DIRECT',RECL=19,FORM='FORMATTED')                    Abre el archivo de datos                
        READ(mtrstr%flind,'(SP,ES18.11E2,A1)',REC=l)vl,enter                                                                        Lee el valor a retornar
        CLOSE(mtrstr%flind)                                                                                                         Cierra el archivo de datos
    ELSE
        CALL error('gt_mtrstrvl_sp: non supported indexes combination')                                                             Sale con error si no se soporta la combinaci�n de �ndices
    END IF   
    !
END SUBROUTINE gt_mtrstrvl_dp
!
!***********************************************************************************************************************************
!************************************** MUESTREAR UN RASTER EN UN CONJUNTO DE PUNTOS ***********************************************
!********************************** @authores: �lvarez-Villa O.D. & Rend�n-�lvarez J.P. ********************************************
!***********************************************************************************************************************************
!
!DESCRIPCI�N:
!Un conjunto de subrutinas �tiles para encontrar el valor de una variable en un conjunto de localizaciones puntuales, a partir de
!una superficie contenida en una variable tipo Metaraster.
!
!OBSERVACI�N IMPORTANTE:
!Cuando uno de los puntos se encuentra por fuera del dominio espacial del metaraster su valor se hace igual a una cantidad (orstval)
!que debe ser especificada por el usuario.
!
!DATOS DE ENTRADA:
!mtrstr     : Variable tipo metaraster                                                                                  {metaraster}
!vblei      : Columna del Metaraster en la que se almacena la varible a muestrear                                           {entero}
!orstval    : Valor a asignar por fuera del dominio del Metaraster                                                   {entero � real}
!XYZ        : Arreglo de coordenadas (x,y,z) de los puntos de muestreo                                                   {real}(:,:)
!
!DATOS DE SALIDA:
!smpldvls    : Arreglo que almacena los valores muesteados en el orden del arreglo de coordenadas               {entero � real}(:)->
!
SUBROUTINE smpl_mtrstr_pnts_i(mtrstr,vblei,orstval,XYZ,smpldvls)
    !
    USE festellustype
    USE util, ONLY: error
    USE tpvrblchckng, ONLY: chckmtrstr
    USE domainprop, ONLY: gt_dmnbnds,gt_vxlZbnds,gt_vxlYbnds,gt_vxlXbnds
    USE raster, ONLY: gt_mtrstrvl
    !
    IMPLICIT NONE
    !
    INTEGER(I4B), INTENT(IN)                         :: vblei,orstval
    INTEGER(I4B), DIMENSION(:), POINTER, INTENT(OUT) :: smpldvls
    REAL(SP), DIMENSION(:,:), INTENT(IN)             :: XYZ
    TYPE(metaraster), INTENT(IN)                     :: mtrstr
    !
    INTEGER(I4B)    :: i,j,k,l,np,col,row,lev,rstval
    REAL(SP)        :: xi,xf,yi,yf,zi,zf,xli,xlf,yli,ylf,zli,zlf
    !
    !Revisi�n de la existencia y pertinencia de los datos de entrada
    CALL chckmtrstr(mtrstr,0,'smpl_mtrstr_pnts_i')                                                                                  El metaraster
    IF(SIZE(XYZ,2)/=3) CALL error('in smpl_mtrstr_pnts_i: invalid coordinates array')                                               Dimensiones del arreglo de coordenadas
    !Muestreo del metaraster en las coordenadas especificadas
    CALL gt_dmnbnds(mtrstr%dom,.true.,xli,yli,zli,xlf,ylf,zlf)                                                                      Obtiene los l�mites espaciales del dominio
    np=SIZE(XYZ,1)                                                                                                                  Calcula el n�mero de puntos a muestrear
    ALLOCATE(smpldvls(np))
    DO l=1,np                                                                                                                       Se busca punto por punto
        IF ((XYZ(l,1)<xli).OR.(XYZ(l,1)>xlf).OR.(XYZ(l,2)<yli).OR.(XYZ(l,2)>ylf).OR.(XYZ(l,3)<zli).OR.(XYZ(l,3)>zlf)) THEN          Si se cumple �sta condici�n el punto se encuentra por fuera del dominio
            smpldvls(l)=orstval                                                                                                     Por lo tanto se asigna el valor que corresponde
        ELSE                                                                                                                        Sino, se b�sca el v�xel
            col=0;row=0;lev=0
            DO k=1,mtrstr%dom%lvls
                CALL gt_vxlZbnds(mtrstr%dom,.true.,k,zi,zf)
                IF ((zi<=XYZ(l,3)).AND.(zf>=XYZ(l,3))) THEN
                    lev=k
                    DO j=1,mtrstr%dom%rws                                                                                                
                        CALL gt_vxlYbnds(mtrstr%dom,.true.,j,yi,yf)
                        IF ((yi<=XYZ(l,2)).AND.(yf>=XYZ(l,2))) THEN
                            row=j
                            DO i=1,mtrstr%dom%clmns
                                CALL gt_vxlXbnds(mtrstr%dom,.true.,i,xi,xf)
                                IF ((xi<=XYZ(l,1)).AND.(xf>=XYZ(l,1))) THEN
                                    col=i
                                    EXIT
                                END IF
                            END DO
                            EXIT
                        END IF
                    END DO
                    EXIT
                END IF
            END DO
            CALL gt_mtrstrvl(mtrstr,.true.,rstval,i,j,k,vblei)                                                                      Se eval�a el Metaraster
            IF (rstval/=mtrstr%falti) THEN                                                                                          
                smpldvls(l)=rstval                                                                                                  Cuando el dato realmente existe
            ELSE
                smpldvls(l)=orstval                                                                                                 En el caso contrario
            END IF
        END IF
    END DO
    !
END SUBROUTINE smpl_mtrstr_pnts_i
!!
SUBROUTINE smpl_mtrstr_pnts_sp(mtrstr,vblei,orstval,XYZ,smpldvls)
    !
    USE festellustype
    USE util, ONLY: error
    USE tpvrblchckng, ONLY: chckmtrstr
    USE domainprop, ONLY: gt_dmnbnds,gt_vxlZbnds,gt_vxlYbnds,gt_vxlXbnds
    USE raster, ONLY: gt_mtrstrvl
    !
    IMPLICIT NONE
    !
    INTEGER(I4B), INTENT(IN)                         :: vblei
    REAL(SP), INTENT(IN)                             :: orstval
    REAL(SP), DIMENSION(:), POINTER, INTENT(OUT)     :: smpldvls
    REAL(SP), DIMENSION(:,:), INTENT(IN)             :: XYZ
    TYPE(metaraster), INTENT(IN)                     :: mtrstr
    !
    INTEGER(I4B)    :: i,j,k,l,np,col,row,lev
    REAL(SP)        :: xi,xf,yi,yf,zi,zf,xli,xlf,yli,ylf,zli,zlf,rstval
    !
    !Revisi�n de la existencia y pertinencia de los datos de entrada
    CALL chckmtrstr(mtrstr,1,'smpl_mtrstr_pnts_sp')                                                                                 El metaraster
    IF(SIZE(XYZ,2)/=3) CALL error('in smpl_mtrstr_pnts_sp: invalid coordinates array')                                              Dimensiones del arreglo de coordenadas
    !Muestreo del metaraster en las coordenadas especificadas
    CALL gt_dmnbnds(mtrstr%dom,.true.,xli,yli,zli,xlf,ylf,zlf)                                                                      Obtiene los l�mites espaciales del dominio
    np=SIZE(XYZ,1)                                                                                                                  Calcula el n�mero de puntos a muestrear
    ALLOCATE(smpldvls(np))
    DO l=1,np                                                                                                                       Se busca punto por punto
        IF ((XYZ(l,1)<xli).OR.(XYZ(l,1)>xlf).OR.(XYZ(l,2)<yli).OR.(XYZ(l,2)>ylf).OR.(XYZ(l,3)<zli).OR.(XYZ(l,3)>zlf)) THEN          Si se cumple �sta condici�n el punto se encuentra por fuera del dominio
            smpldvls(l)=orstval                                                                                                     Por lo tanto se asigna el valor que corresponde
        ELSE                                                                                                                        Sino, se b�sca el v�xel
            col=0;row=0;lev=0
            DO k=1,mtrstr%dom%lvls
                CALL gt_vxlZbnds(mtrstr%dom,.true.,k,zi,zf)
                IF ((zi<=XYZ(l,3)).AND.(zf>=XYZ(l,3))) THEN
                    lev=k
                    DO j=1,mtrstr%dom%rws                                                                                                
                        CALL gt_vxlYbnds(mtrstr%dom,.true.,j,yi,yf)
                        IF ((yi<=XYZ(l,2)).AND.(yf>=XYZ(l,2))) THEN
                            row=j
                            DO i=1,mtrstr%dom%clmns
                                CALL gt_vxlXbnds(mtrstr%dom,.true.,i,xi,xf)
                                IF ((xi<=XYZ(l,1)).AND.(xf>=XYZ(l,1))) THEN
                                    col=i
                                    EXIT
                                END IF
                            END DO
                            EXIT
                        END IF
                    END DO
                    EXIT
                END IF
            END DO
            CALL gt_mtrstrvl(mtrstr,.true.,rstval,i,j,k,vblei)                                                                      Se eval�a el Metaraster
            IF (rstval/=mtrstr%falts) THEN                                                                                           
                smpldvls(l)=rstval                                                                                                  Cuando el dato realmente existe
            ELSE
                smpldvls(l)=orstval                                                                                                 En el caso contrario
            END IF
        END IF
    END DO
    !
END SUBROUTINE smpl_mtrstr_pnts_sp
!!
SUBROUTINE smpl_mtrstr_pnts_dp(mtrstr,vblei,orstval,XYZ,smpldvls)
    !
    USE festellustype
    USE util, ONLY: error
    USE tpvrblchckng, ONLY: chckmtrstr
    USE domainprop, ONLY: gt_dmnbnds,gt_vxlZbnds,gt_vxlYbnds,gt_vxlXbnds
    USE raster, ONLY: gt_mtrstrvl
    !
    IMPLICIT NONE
    !
    INTEGER(I4B), INTENT(IN)                         :: vblei
    REAL(DP), INTENT(IN)                             :: orstval
    REAL(DP), DIMENSION(:), POINTER, INTENT(OUT)     :: smpldvls
    REAL(DP), DIMENSION(:,:), INTENT(IN)             :: XYZ
    TYPE(metaraster), INTENT(IN)                     :: mtrstr
    !
    INTEGER(I4B)    :: i,j,k,l,np,col,row,lev
    REAL(DP)        :: xi,xf,yi,yf,zi,zf,xli,xlf,yli,ylf,zli,zlf,rstval
    !
    !Revisi�n de la existencia y pertinencia de los datos de entrada
    CALL chckmtrstr(mtrstr,2,'smpl_mtrstr_pnts_dp')                                                                                 El metaraster
    IF(SIZE(XYZ,2)/=3) CALL error('in smpl_mtrstr_pnts_dp: invalid coordinates array')                                              Dimensiones del arreglo de coordenadas
    !Muestreo del metaraster en las coordenadas especificadas
    CALL gt_dmnbnds(mtrstr%dom,.true.,xli,yli,zli,xlf,ylf,zlf)                                                                      Obtiene los l�mites espaciales del dominio
    np=SIZE(XYZ,1)                                                                                                                  Calcula el n�mero de puntos a muestrear
    ALLOCATE(smpldvls(np))
    DO l=1,np                                                                                                                       Se busca punto por punto
        IF ((XYZ(l,1)<xli).OR.(XYZ(l,1)>xlf).OR.(XYZ(l,2)<yli).OR.(XYZ(l,2)>ylf).OR.(XYZ(l,3)<zli).OR.(XYZ(l,3)>zlf)) THEN          Si se cumple �sta condici�n el punto se encuentra por fuera del dominio
            smpldvls(l)=orstval                                                                                                     Por lo tanto se asigna el valor que corresponde
        ELSE                                                                                                                        Sino, se b�sca el v�xel
            col=0;row=0;lev=0
            DO k=1,mtrstr%dom%lvls
                CALL gt_vxlZbnds(mtrstr%dom,.true.,k,zi,zf)
                IF ((zi<=XYZ(l,3)).AND.(zf>=XYZ(l,3))) THEN
                    lev=k
                    DO j=1,mtrstr%dom%rws                                                                                                
                        CALL gt_vxlYbnds(mtrstr%dom,.true.,j,yi,yf)
                        IF ((yi<=XYZ(l,2)).AND.(yf>=XYZ(l,2))) THEN
                            row=j
                            DO i=1,mtrstr%dom%clmns
                                CALL gt_vxlXbnds(mtrstr%dom,.true.,i,xi,xf)
                                IF ((xi<=XYZ(l,1)).AND.(xf>=XYZ(l,1))) THEN
                                    col=i
                                    EXIT
                                END IF
                            END DO
                            EXIT
                        END IF
                    END DO
                    EXIT
                END IF
            END DO
            CALL gt_mtrstrvl(mtrstr,.true.,rstval,i,j,k,vblei)                                                                      Se eval�a el Metaraster
            IF (rstval/=mtrstr%faltd) THEN                                                                                           
                smpldvls(l)=rstval                                                                                                  Cuando el dato realmente existe
            ELSE
                smpldvls(l)=orstval                                                                                                 En el caso contrario
            END IF
        END IF
    END DO
    !
END SUBROUTINE smpl_mtrstr_pnts_dp