ALL: ANISOFLOW
   CFLAGS     =
   FFLAGS     = 
   CPPFLAGS   =
   FPPFLAGS   =
   CLEANFILES = ANISOFLOW

ANISOFLOW: 01_Types.o 02_Interface.o 03_View.o 04_Geometry.o 05_Properties.o 06_BoundaryConditions.o 07_BuildSystem.o 08_Solver.o ANISOFLOW.o chkopts
	-${FLINKER} ANISOFLOW.o 01_Types.o 02_Interface.o 03_View.o 04_Geometry.o 05_Properties.o 06_BoundaryConditions.o 07_BuildSystem.o 08_Solver.o -o ANISOFLOW ${PETSC_SYS_LIB}
	${RM} -f *.o *mod


include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
