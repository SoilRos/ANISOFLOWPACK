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
MODULE anisotype
    !
    !*******************************************************************************************************************************
    !*******TIPOS DE VARIABLES CONSIDERADAS EN LA GENERACIÓN DE CÓDIGO (Extraida de FDPACK)*****************************************
    !*******@authores: Press et al. (1986)******************************************************************************************
    !*******************************************************************************************************************************
    !
	INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
	INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
	INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
	INTEGER, PARAMETER :: SP = KIND(1.0)
	INTEGER, PARAMETER :: DP = KIND(1.0D0)
	INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
	INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
	INTEGER, PARAMETER :: LGT = KIND(.true.)
    !
    !*******************************************************************************************************************************
    !*************CONSTANTES ÚTILES (Extraida de FDPACK)****************************************************************************
    !*******@authores: Press et al. (1986)******************************************************************************************
    !*******************************************************************************************************************************
	!
	REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
	REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
	REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
	REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
	REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
	REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
	REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
	REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp
	REAL(DP), PARAMETER :: SQRT2_D=1.41421356237309504880168872420969807856967_dp
	REAL(DP), PARAMETER :: EULER_D=0.5772156649015328606065120900824024310422_dp
    !
    !*******************************************************************************************************************************
    !*********VARIABLE TIPO PARA LAS MATRICES DISPERSAS EN COORDENADAS**************************************************************
    !*******@authores: Press et al. (1986) @modified: Alvarez-Villa O.D.************************************************************
    !*******************************************************************************************************************************
    !    
    !n=Longitud de la amtriz dispersa 
    !val=Valor del coeficiencte en la matriz dispersa
    !irow=Coordenada en x de la matriz densa para la dispersa
    !jcol=Coordenada en y de la matriz densa para la dispersa
    !
    TYPE dispersed
        INTEGER(I4B)                        :: n,len
	REAL(DP), DIMENSION(:), POINTER     :: val
	INTEGER(I4B), DIMENSION(:), POINTER :: irow
	INTEGER(I4B), DIMENSION(:), POINTER :: jcol
    END TYPE dispersed
    !    
    !*******************************************************************************************************************************
    !**********VARIBALE TIPO PARA LA GEOMETRIA DEL ACUIFERO*************************************************************************
    !**********@authores: Alvarez-Villa O.D., Rendon-Avarez J.P.********************************************************************
    !*******************************************************************************************************************************
    !
    !Variables tipo para la geometria del dominio espacial. Los parametros asociados son:
    !exist          :Variable logica que indica si el dominio se ha construido
    !infile         :Variable logica que indica si el dataval esta en archivo o en memoria (.false.=RAM  .true.=disco)
    !flind          :Entero que indica la localizacion en memoria del archivo de datos si infile=true
    !X              :Coordenadas de los bordes de las columnas de la malla s:flotante, d:doble
    !Y              :Coordenadas de los bordes de las filas de la malla  d:doble
    !Z              :Coordenadas de los bordes de las capas de la malla  d:doble
    !rws            :Numero de filas en los cuales se ha dividido el dominio
    !clmns          :Numero de columnas en las cuales se ha dividido el dominio
    !lvls           :Numero de niveles en los cuales se ha dividido el dominio  
    !
    TYPE domain
        LOGICAL(LGT)                    :: exist=.false.,infile=.false.
        INTEGER(I4B)                    :: rws,clmns,lvls,flind=-9999 
        CHARACTER(200)                  :: rtdom
        REAL(SP), DIMENSION(:), POINTER :: Xs,Ys,Zs
        REAL(DP), DIMENSION(:), POINTER :: Xd,Yd,Zd
    END TYPE domain

    !
    !*******************************************************************************************************************************
    !**************VARIABLE TIPO PARA LA TOPOLOGIA DEL ACUIFERO*********************************************************************
    !**************@autores: Alvarez-Villa O.D., Perez Kevin************************************************************************
    !*******************************************************************************************************************************
    !
    !Variable tipo que contiene la topologia del acuifero. Los parametros asociados son:
    !act: Numero de celdas activas.
    !cth: Numero de celdas de nivel impuesto.
    !dph: Numero de celdas de nivel externo.
    !hastopo: Variable logica que indica existencia de variable tipo.
    !ijactv: Arreglo con numeracion topologica de las celdas activas efectivas.
    !ijcth: Arreglo con numeracion topologica de las celdas de condicion de contorno de altura piezometrica constante. 
    !ijdph: Arreglo con numeracion topologica de las celdas con condicion de contorno dependientes de una altura externa.
    !
    TYPE topology
        INTEGER(I4B)                   :: act,cth,dph
        LOGICAL(LGT)                   :: hastopo=.false.
        INTEGER,DIMENSION(:,:),POINTER :: ijactv,ijcth,ijdph
    END TYPE topology  
    
END MODULE anisotype