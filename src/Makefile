Run.exe: 01_Types.o 02_Util.o 03_Type_Checking.o 04_Geometry.o 05_IOdata.o 06_K_operations.o 07_Topology.o 08_IsoSystem.o 09_AnisoSystem.o 11_Simulation.o 12_LinearSolvers.o Run.o chkopts
	-${FLINKER} Run.o 01_Types.o 02_Util.o 03_Type_Checking.o 04_Geometry.o 05_IOdata.o 06_K_operations.o 07_Topology.o 08_IsoSystem.o 09_AnisoSystem.o 11_Simulation.o 12_LinearSolvers.o -o Run.exe ${PETSC_SYS_LIB}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

#Run.o: Run.f90
#	gfortran -c Run.f90

01_Types.o: 01_Types.f90
	gfortran -c 01_Types.f90

02_Util.o: 02_Util.f90
	gfortran -c 02_Util.f90

03_Type_Checking.o: 03_Type_Checking.f90
	gfortran -c 03_Type_Checking.f90

04_Geometry.o: 04_Geometry.f90
	gfortran -c 04_Geometry.f90

05_IOdata.o: 05_IOdata.f90
	gfortran -c 05_IOdata.f90

06_K_operations.o: 06_K_operations.f90
	gfortran -c 06_K_operations.f90

07_Topology.o: 07_Topology.f90
	gfortran -c 07_Topology.f90

08_IsoSystem.o: 08_IsoSystem.f90
	gfortran -c 08_IsoSystem.f90

09_AnisoSystem.o: 09_AnisoSystem.f90
	gfortran -c 09_AnisoSystem.f90

11_Simulation.o: 11_Simulation.f90
	gfortran -c 11_Simulation.f90

12_LinearSolvers.o:	12_LinearSolvers.f90
	gfortran -c 12_LinearSolvers.f90
