ANISOFLOW.exe: 01_Types.o 02_Interface.o 03_Geometry.o 04_Properties.o 05_BoundaryConditions.o 06_BuildSystem.o 07_Solver.o ANISOFLOW.o chkopts
	-${FLINKER} ANISOFLOW.o 01_Types.o 02_Interface.o 03_Geometry.o 04_Properties.o 05_BoundaryConditions.o 06_BuildSystem.o 07_Solver.o -o ANISOFLOW.exe ${PETSC_SYS_LIB}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules