PROGRAM ANISOFLOW_RANDOM_FIELD_GENERATOR

    IMPLICIT NONE

    REAL                    :: RealWidth(3),IsoSizeReal,OAnisoSizeReal,NOAnisoSizeReal
    INTEGER                 :: IndexWidth(3),IsoSize,OAnisoSize,NOAnisoSize,ZonesSize,SimulationDays
    INTEGER                 :: NArgs,i
    CHARACTER(LEN=20)       :: CharIndexWidth(3)
    CHARACTER(LEN=200)      :: ProjectName,Arg

    RealWidth(:)=400.D0
    IndexWidth(:)=25
    SimulationDays=365*4

    NArgs=COMMAND_ARGUMENT_COUNT()
    IF (Nargs.GT.0) THEN
        DO i=1,NArgs
            CALL GET_COMMAND_ARGUMENT(i,Arg)
            READ(Arg,"(I5)") IndexWidth(i)
            IF (i.EQ.3) EXIT
        END DO
    END IF

    CALL RANDOM_NUMBER(IsoSizeReal);        IsoSize=INT(3*IsoSizeReal)
    CALL RANDOM_NUMBER(OAnisoSizeReal);     OAnisoSize=INT(0*OAnisoSizeReal)
    CALL RANDOM_NUMBER(NOAnisoSizeReal);    NOAnisoSize=INT(2*NOAnisoSizeReal)
    ZonesSize=IsoSize+OAnisoSize+NOAnisoSize

    WRITE(CharIndexWidth(1),*)IndexWidth(1)
    WRITE(CharIndexWidth(2),*)IndexWidth(2)
    WRITE(CharIndexWidth(3),*)IndexWidth(3)
    ProjectName=ADJUSTL(TRIM(CharIndexWidth(1))//"x")
    ProjectName=TRIM(ProjectName)//ADJUSTL(TRIM(CharIndexWidth(2))//"x")
    ProjectName=TRIM(ProjectName)//ADJUSTL(TRIM(CharIndexWidth(3)))

    CALL CreateGmtry(ProjectName,RealWidth,IndexWidth)
    CALL CreateProperties(ProjectName,IsoSize,OAnisoSize,NOAnisoSize)
    CALL CreateZones(ProjectName,IndexWidth,ZonesSize)
    CALL CreateBC(ProjectName,IndexWidth,SimulationDays)
    CALL CreateExecutionLine(ProjectName)
    CALL CreatePerformanceTest(ProjectName)

END PROGRAM ANISOFLOW_RANDOM_FIELD_GENERATOR

SUBROUTINE CreateGmtry(ProjectName,RealWidth,IndexWidth)

    IMPLICIT NONE

    REAL,INTENT(IN)                 :: RealWidth(3)
    INTEGER,INTENT(IN)              :: IndexWidth(3)
    CHARACTER(LEN=200),INTENT(IN)   :: ProjectName


    CHARACTER(LEN=200)              :: FileName
    INTEGER                         :: i
    INTEGER,PARAMETER               :: u=01

    FileName=TRIM(ProjectName)//"_gmtry.txt"
    OPEN(u,FILE=FileName, STATUS="REPLACE")
    
    WRITE(u,*)"COLUMNS  ",IndexWidth(1)
    WRITE(u,*)"ROWS     ",IndexWidth(2)
    WRITE(u,*)"LEVELS   ",IndexWidth(3)

    DO i=0,IndexWidth(1)
        WRITE(u,*)i*RealWidth(1)/IndexWidth(1)
    END DO
    DO i=0,IndexWidth(2)
        WRITE(u,*)i*RealWidth(2)/IndexWidth(2)
    END DO
    DO i=0,IndexWidth(3)
        WRITE(u,*)i*RealWidth(3)/IndexWidth(3)
    END DO
    CLOSE(u)


END SUBROUTINE CreateGmtry

SUBROUTINE CreateProperties(ProjectName,IsoSize,OAnisoSize,NOAnisoSize)

    IMPLICIT NONE

    INTEGER,INTENT(IN)              :: IsoSize,OAnisoSize,NOAnisoSize
    CHARACTER(LEN=200),INTENT(IN)   :: ProjectName

    CHARACTER(LEN=200)              :: FileName
    REAL                            :: RandCvt(7)
    INTEGER                         :: i
    INTEGER,PARAMETER               :: u=01

    FileName=TRIM(ProjectName)//"_properties.txt"
    OPEN(u,FILE=FileName, STATUS="REPLACE")

    DO i=1,(IsoSize+OAnisoSize+NOAnisoSize)
        CALL RANDOM_NUMBER(RandCvt)
        WRITE(u,*)"!------------------------------------------"
        WRITE(u,*)"Material: ",i
        WRITE(u,*)" "

        IF (i.LE.IsoSize) THEN
            WRITE(u,*)"k isotropic"
            WRITE(u,*) 5*RandCvt(1)
        ELSEIF (i.LE.(IsoSize+OAnisoSize)) THEN
            WRITE(u,*)"k anisotropic"
            WRITE(u,*) 5*RandCvt(1),5*RandCvt(2),5*RandCvt(3)
        ELSE
            WRITE(u,*)"k non-orthogonal anisotropic"
            WRITE(u,*)5*RandCvt(1),5*RandCvt(2),5*RandCvt(3),2*RandCvt(4),0.1*RandCvt(5),0.1*RandCvt(6)
        END IF

        WRITE(u,*)" "
        WRITE(u,*)"specific storage"
        WRITE(u,*)RandCvt(7)    !!!!!!!!!!!
        WRITE(u,*)" "
        WRITE(u,*)"longitudinal dispersivity"
        WRITE(u,*)"10.0"
        WRITE(u,*)" "
        WRITE(u,*)"transverse dispersivity"
        WRITE(u,*)"1.0d0"
        WRITE(u,*)" "
        WRITE(u,*)"vertical transverse dispersivity"
        WRITE(u,*)"1.0d0"
        WRITE(u,*)" "
        WRITE(u,*)"porosity"
        WRITE(u,*)"0.005d0"
        WRITE(u,*)" "
        WRITE(u,*)"tortuosity"
        WRITE(u,*)"0.1"
        WRITE(u,*)" "
        WRITE(u,*)"bulk density"
        WRITE(u,*)"2650.0d0"
        WRITE(u,*)" "
        WRITE(u,*)"end material"
        WRITE(u,*)" "
    END DO
    
    CLOSE(u)

END SUBROUTINE CreateProperties

SUBROUTINE CreateZones(ProjectName,IndexWidth,ZonesSize)

    IMPLICIT NONE

    INTEGER,INTENT(IN)              :: ZonesSize,IndexWidth(3)
    CHARACTER(LEN=200),INTENT(IN)   :: ProjectName

    CHARACTER(LEN=200)              :: FileName
    REAL                            :: ZoneReal,RandGroup,RandGroupX,RandGroupY,RandGroupZ
    INTEGER                         :: i,j,k,GroupX,GroupY,GroupZ
    INTEGER,ALLOCATABLE             :: ZoneInt(:,:,:)
    INTEGER,PARAMETER               :: u=01


    ALLOCATE(ZoneInt(IndexWidth(1),IndexWidth(2),IndexWidth(3)))

    FileName=TRIM(ProjectName)//"_zones.txt"
    OPEN(u,FILE=FileName, STATUS="REPLACE")

    CALL RANDOM_NUMBER(ZoneReal);ZoneInt=INT((ZonesSize)*ZoneReal+1)
    ZoneInt(:,:,:)=ZoneInt

    DO i=1,IndexWidth(1)
        DO j=1,IndexWidth(2)
            DO k=1,IndexWidth(3)
                CALL RANDOM_NUMBER(ZoneReal)
                CALL RANDOM_NUMBER(RandGroup)
                RandGroupX=0.D0;RandGroupY=0.D0;RandGroupZ=0.D0
                IF (i.GT.1) CALL RANDOM_NUMBER(RandGroupX)
                IF (j.GT.1) CALL RANDOM_NUMBER(RandGroupY)
                IF (k.GT.1) CALL RANDOM_NUMBER(RandGroupZ)
                GroupX=i-INT(2*RandGroupX)
                GroupY=j-INT(2*RandGroupY)
                GroupZ=k-INT(2*RandGroupZ)
                ZoneInt(i,j,k)=INT((ZonesSize)*ZoneReal+1)
                ZoneInt(i,j,k)=ZoneInt(GroupX,GroupY,GroupZ)
                WRITE(u,*)ZoneInt(i,j,k)
            END DO
        END DO
    END DO

END SUBROUTINE CreateZones

SUBROUTINE CreateBC(ProjectName,IndexWidth,SimulationDays)

    IMPLICIT NONE

    INTEGER,INTENT(IN)              :: IndexWidth(3),SimulationDays
    CHARACTER(LEN=200),INTENT(IN)   :: ProjectName

    CHARACTER(LEN=200)              :: FileName
    INTEGER                         :: i,j,k,DirichSize
    LOGICAL                         :: StressSinkPeriod,StressSourcePeriod

    INTEGER,PARAMETER               :: u=01
    REAL,PARAMETER                  :: PI = 3.1415927


    FileName=TRIM(ProjectName)//"_bc.txt"
    OPEN(u,FILE=FileName, STATUS="REPLACE")

    WRITE(u,*)"Timezones: ",SimulationDays

    DirichSize=2*IndexWidth(1)*IndexWidth(3)
    i=1
    WRITE(u,*)"Timezone(",i,"):",1
    WRITE(u,*)"DT:",1.0
    WRITE(u,*)"DirichletBoundaryCondition(",i,"):", DirichSize
    DO j=1,IndexWidth(1)
        DO k=1,IndexWidth(3)
            WRITE(u,*)j,1            ,k,10*(COS(2*PI*(j-1)/(IndexWidth(1)-1)+PI)*0.25+0.75)*(COS(2*PI*(i-1)/(365-1))*0.5+0.5)
            WRITE(u,*)j,IndexWidth(2),k,10*(COS(2*PI*(j-1)/(IndexWidth(1)-1)+PI)*0.25+0.75)*(COS(2*PI*(i-1)/(365-1))*0.5+0.5)
        END DO
    END DO
    WRITE(u,*)"SourceBoundaryCondition(",i,"):",0
    WRITE(u,*)"CauchyBoundaryCondition(",i,"):",0

    DO i=2,SimulationDays
        WRITE(u,*)"Timezone(",i,"):",1
        WRITE(u,*)"DT:",1.0

        StressSourcePeriod=(((i/(2*365.D0))-FLOOR(i/(2*365.D0))).GT.0.5D0).AND.(((i/(2*365.D0))-FLOOR(i/(2*365.D0))).LT.0.75D0)
        StressSinkPeriod=  (((i/(2*365.D0))-FLOOR(i/(2*365.D0))).LT.0.25D0)
        DirichSize=2*IndexWidth(1)*IndexWidth(3)
        IF (StressSinkPeriod.OR.StressSourcePeriod) DirichSize=DirichSize+3
        WRITE(u,*)"DirichletBoundaryCondition(",i,"):", DirichSize
        DO j=1,IndexWidth(1)
            DO k=1,IndexWidth(3)
                WRITE(u,*)j,1            ,k,10*(COS(2*PI*(j-1)/(IndexWidth(1)-1)+PI)*0.25+0.75)*(COS(2*PI*(i-1)/(365-1))*0.5+0.5)
                WRITE(u,*)j,IndexWidth(2),k,10*(COS(2*PI*(j-1)/(IndexWidth(1)-1)+PI)*0.25+0.75)*(COS(2*PI*(i-1)/(365-1))*0.5+0.5)
            END DO
        END DO
        IF (StressSinkPeriod) THEN
            WRITE(u,*)INT(IndexWidth(1)/2),INT(IndexWidth(2)/2),1,0
            WRITE(u,*)INT(IndexWidth(1)/2),INT(IndexWidth(2)/2),INT(IndexWidth(3)/2),0
            WRITE(u,*)INT(IndexWidth(1)/2),INT(IndexWidth(2)/2),INT(IndexWidth(3)),0
        ELSEIF (StressSourcePeriod) THEN
            WRITE(u,*)INT(IndexWidth(1)/2),INT(IndexWidth(2)/2),1,20
            WRITE(u,*)INT(IndexWidth(1)/2),INT(IndexWidth(2)/2),INT(IndexWidth(3)/2),20
            WRITE(u,*)INT(IndexWidth(1)/2),INT(IndexWidth(2)/2),INT(IndexWidth(3)),20
        END IF
        WRITE(u,*)"SourceBoundaryCondition(",i,"):",0
        WRITE(u,*)"CauchyBoundaryCondition(",i,"):",0
    END DO

    CLOSE(u)

END SUBROUTINE CreateBC

SUBROUTINE CreateExecutionLine(ProjectName)

    IMPLICIT NONE

    CHARACTER(LEN=200),INTENT(IN)   :: ProjectName

    CHARACTER(LEN=200)              :: FileName
    INTEGER,PARAMETER               :: u=01

    FileName=TRIM(ProjectName)//".sh"
    OPEN(u,FILE=FileName, STATUS="REPLACE")

    WRITE(u,*)"../../src/ANISOFLOW \"
    WRITE(u,*)"-Input_type_gmtry 2 \"
    WRITE(u,*)"-Input_file_gmtry ",TRIM(ProjectName)//"_gmtry.txt \"
    WRITE(u,*)"-Input_type_tplgy 0 \"
    WRITE(u,*)"-Input_file_ppt ",TRIM(ProjectName)//"_properties.txt \"
    WRITE(u,*)"-Input_file_ppt_by_zones ",TRIM(ProjectName)//"_zones.txt \"
    WRITE(u,*)"-Input_file_bc ",TRIM(ProjectName)//"_bc.txt \"
    WRITE(u,*)"-Input_type_bc 2 \"
    WRITE(u,*)"-Run_options_scheme 2 \"
    WRITE(u,*)"-Output_type_ppt 3 \"
    WRITE(u,*)"-Project_name ",TRIM(ProjectName)," "
    WRITE(u,*)" "

    WRITE(u,*)"../../src/ANISOFLOW \"
    WRITE(u,*)"-Input_type_gmtry 2 \"
    WRITE(u,*)"-Input_file_gmtry ",TRIM(ProjectName)//"_gmtry.txt \"
    WRITE(u,*)"-Input_type_tplgy 0 \"
    WRITE(u,*)"-Input_file_ppt ",TRIM(ProjectName)//"_properties.txt \"
    WRITE(u,*)"-Input_file_ppt_by_zones ",TRIM(ProjectName)//"_zones.txt \"
    WRITE(u,*)"-Input_file_bc ",TRIM(ProjectName)//"_bc.txt \"
    WRITE(u,*)"-Input_type_bc 2 \"
    WRITE(u,*)"-Run_options_scheme 2 \"
    WRITE(u,*)"-Run_options_time 1 \"
    WRITE(u,*)"-Output_type_ppt 3 \"
    WRITE(u,*)"-Project_name ",TRIM(ProjectName)," \"
    WRITE(u,*)"-Input_type_init_sol 3 \"
    WRITE(u,*)"-Input_file_init_sol ",TRIM(ProjectName)//"_Sol.h5 "

    CLOSE(u)
    CALL CHMOD(FileName,"755")

END SUBROUTINE CreateExecutionLine

SUBROUTINE CreatePerformanceTest(ProjectName)

    IMPLICIT NONE

    CHARACTER(LEN=200),INTENT(IN)   :: ProjectName

    CHARACTER(LEN=200)              :: FileName
    INTEGER,PARAMETER               :: u=01

    FileName=TRIM(ProjectName)//"_PerformanceTest.sh"
    OPEN(u,FILE=FileName, STATUS="REPLACE")

    WRITE(u,*)"# Execution Line must be: ./"//TRIM(FileName)//" np_max"
    WRITE(u,*)"np_max= $1"
    WRITE(u,*)"declare -a KSP=(""cg"" ""gmres"")"
    WRITE(u,*)"declare -a PC=(""none"" ""jacobi"" ""bjacobi"")"

    WRITE(u,*)"for np in `seq 1 ${np_max}`; do"
    WRITE(u,*)"    for ksptype in ""${KSP[@]}""; do"
    WRITE(u,*)"        for pctype in ""${PC[@]}""; do"
    WRITE(u,*)"            echo ""*************************************************"""
    WRITE(u,*)"            echo ""*************** Beginning new run ***************"""
    WRITE(u,*)"            echo ""Options"""
    WRITE(u,*)"            echo ""np: $np"""
    WRITE(u,*)"            echo ""Solver: $ksptype"""
    WRITE(u,*)"            echo ""Preconditioner: $pctype"""
    WRITE(u,*)"            echo ""************** Initializing program *************"""
    WRITE(u,*)"            mpiexec -n ${np} \"
    WRITE(u,*)"            ../../src/ANISOFLOW \"
    WRITE(u,*)"            -Input_type_gmtry 2 \"
    WRITE(u,*)"            -Input_file_gmtry ",TRIM(ProjectName)//"_gmtry.txt \"
    WRITE(u,*)"            -Input_type_tplgy 0 \"
    WRITE(u,*)"            -Input_file_ppt ",TRIM(ProjectName)//"_properties.txt \"
    WRITE(u,*)"            -Input_file_ppt_by_zones ",TRIM(ProjectName)//"_zones.txt \"
    WRITE(u,*)"            -Input_file_bc ",TRIM(ProjectName)//"_bc.txt \"
    WRITE(u,*)"            -Input_type_bc 2 \"
    WRITE(u,*)"            -Run_options_scheme 2 \"
    WRITE(u,*)"            -Output_type_ppt 3 \"
    WRITE(u,*)"            -Project_name ",TRIM(ProjectName),"\"
    WRITE(u,*)"            -ksp_type ${ksptype} \"
    WRITE(u,*)"            -pc_type  ${pctype} \"
    WRITE(u,*)"             "

    WRITE(u,*)"            mpiexec -n ${np} \"
    WRITE(u,*)"            ../../src/ANISOFLOW \"
    WRITE(u,*)"            -Input_type_gmtry 2 \"
    WRITE(u,*)"            -Input_file_gmtry ",TRIM(ProjectName)//"_gmtry.txt \"
    WRITE(u,*)"            -Input_type_tplgy 0 \"
    WRITE(u,*)"            -Input_file_ppt ",TRIM(ProjectName)//"_properties.txt \"
    WRITE(u,*)"            -Input_file_ppt_by_zones ",TRIM(ProjectName)//"_zones.txt \"
    WRITE(u,*)"            -Input_file_bc ",TRIM(ProjectName)//"_bc.txt \"
    WRITE(u,*)"            -Input_type_bc 2 \"
    WRITE(u,*)"            -Run_options_scheme 2 \"
    WRITE(u,*)"            -Run_options_time 1 \"
    WRITE(u,*)"            -Output_type_ppt 3 \"
    WRITE(u,*)"            -Project_name ",TRIM(ProjectName)," \"
    WRITE(u,*)"            -Input_type_init_sol 3 \"
    WRITE(u,*)"            -Input_file_init_sol ",TRIM(ProjectName)//"_Sol.h5 \"
    WRITE(u,*)"            -ksp_type ${ksptype} \"
    WRITE(u,*)"            -pc_type  ${pctype} \"
    WRITE(u,*)"            -log_view ascii:"//TRIM(ProjectName)//"_NP${np}_KSP${ksptype}_PC${pctype}.log \"
    WRITE(u,*)"            | tee "//TRIM(ProjectName)//"_NP${np}_KSP${ksptype}_PC${pctype}.log_term "
    WRITE(u,*)"            echo ""***************** Ending program ****************"""
    WRITE(u,*)"        done"
    WRITE(u,*)"    done"
    WRITE(u,*)"done"
    CLOSE(u)
    CALL CHMOD(FileName,"755")

END SUBROUTINE CreatePerformanceTest