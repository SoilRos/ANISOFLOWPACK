ALL: ANISOFLOW
   CFLAGS     =
   FFLAGS     = 
   CPPFLAGS   =
   FPPFLAGS   =
   CLEANFILES = ANISOFLOW

ANISOFLOW_OBJ = 01_Types.o 02_Operators.o 03_Interface.o 04_View.o 	 \
				05_Geometry.o 06_Properties.o 07_BoundaryConditions.o\
				08_BuildSystem.o 09_Solver.o ANISOFLOW.o

ANISOFLOW: ${ANISOFLOW_OBJ} chkopts
	-${FLINKER} ${ANISOFLOW_OBJ} -o ANISOFLOW ${PETSC_SYS_LIB}
	${RM} -f *.o *mod


include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
