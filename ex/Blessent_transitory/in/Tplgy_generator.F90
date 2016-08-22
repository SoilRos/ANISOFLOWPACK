PROGRAM BC_GENERATOR
    IMPLICIT NONE
    INTEGER                 :: i,j,k,l,u,sizeBC,cont,ValI
    INTEGER,ALLOCATABLE     :: BC_index(:)
    REAL,ALLOCATABLE        :: BC_value(:)

    u=100

    sizeBC=320
    ALLOCATE(BC_index(sizeBC))
    ALLOCATE(BC_value(sizeBC))
    OPEN(u,FILE='BC_lateralboundary.txt', STATUS='OLD')
    DO i=1,sizeBC
        READ(u,*) BC_index(i),BC_value(i)
    END DO
    CLOSE(u)

    OPEN(u,FILE='ANISOFLOW_TEST.tplgy', STATUS='NEW')
    cont=1

    DO k=1,41
        DO j=1,81
            DO i=1,81
                ValI=1
                DO l=1,sizeBC
                    IF (BC_index(l).EQ.cont) THEN
                        ValI=2
                        EXIT
                    END IF
                END DO
                WRITE(u,*) ValI
                cont=cont+1
            END DO
        END DO
    END DO
    CLOSE(u)

END PROGRAM BC_GENERATOR