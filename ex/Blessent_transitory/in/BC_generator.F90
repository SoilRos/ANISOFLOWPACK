PROGRAM BC_GENERATOR
    IMPLICIT NONE
    INTEGER                 :: i,j,t,u,dtsize,timezones,sizeBC,x,y,z
    INTEGER,ALLOCATABLE     :: BC_index(:)
    REAL,ALLOCATABLE        :: BC_value(:),BC_value_tmp(:)
    INTEGER                 :: dt
    REAL, PARAMETER         :: Pi = 3.1415927

    u=100
    i=1
    timezones=50
    dtsize=1
    dt=1

    sizeBC=320
    ALLOCATE(BC_index(sizeBC))
    ALLOCATE(BC_value(sizeBC))
    OPEN(u,FILE='BC_lateralboundary.txt', STATUS='OLD')
    DO i=1,sizeBC
        READ(u,*) BC_index(i),BC_value(i)
    END DO
    CLOSE(u)
    ALLOCATE(BC_value_tmp(sizeBC))
    OPEN(u,FILE='ANISOFLOW_TEST.bc', STATUS='REPLACE')
    WRITE(u,*)'Timezones: ',timezones
    i=1
    DO t=1,dt*timezones,dt
        BC_value_tmp(:)=BC_value(:)*(COS(2*PI*t/365)*0.5+0.5)
        WRITE(u,*)'Timezone( ',i,'): ', dtsize
        WRITE(u,*)'DT: ',dt
        WRITE(u,*)'DirichletBoundaryCondition( ',i,' ): ', sizeBC
        DO j=1,sizeBC
            WRITE(u,*)BC_index(j),BC_value_tmp(j)
        END DO
        WRITE(u,*)'SourceBoundaryCondition( ',i,' ): 4'
        x=20;y=20;z=39
        WRITE(u,*)81*81*(z-1)+81*(y-1)+x,-40 !Esquina 20,20,32 zona 4
        x=61;y=20;z=32
        WRITE(u,*)81*81*(z-1)+81*(y-1)+x,-40 !Esquina 20,48,32 zona 2
        x=61;y=61;z=38
        WRITE(u,*)81*81*(z-1)+81*(y-1)+x,-40 !Esquina 61,62,31 zona 3
        x=20;y=61;z=35
        WRITE(u,*)81*81*(z-1)+81*(y-1)+x,-40 !Esquina 61,62,31 zona 1
        WRITE(u,*)'CauchyBoundaryCondition( ',i,' ): 0'
        i=i+1
    END DO
    CLOSE(u)

END PROGRAM BC_GENERATOR
