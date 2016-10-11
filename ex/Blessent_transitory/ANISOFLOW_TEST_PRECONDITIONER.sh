#!/bin/bash

np_max=4

#PBS -N ANISOFLOW			# nombre identificador del trabajo
#PBS -q cola$np_max			  	 # Nombre de la cola de ejecucion
#PBS -l nodes=1:ppn=$np_max           # nodes=1 numero de nodos a emplear
     				  # ppn=8 numero de procesadores a emplear
#PBS  -l walltime=10000:00:00   # tiempo maximo de ejecucion
#PBS -o salida.OUT              # archivo de salida
#PBS -e salida.ERR              # archivo de error

ProgramDir=../../src/
ProgramName=ANISOFLOW

InputType=1													# Input files type
InputTypeTplgy=1 											# Input toplogy file type
InputTypeGmtry=2 											# Input gemoetry file type
InputTypeBC=1 												# Input boundary conditions file type
InputTypeInitSol=3
InputDir=in/ 												# Input directory
InputFileGmtry=ANISOFLOW_TEST.grid 								# Geometry file
InputFileTplgy=ANISOFLOW_TEST.tplgy 								# Topology file
InputFileCvt=ANISOFLOW_TEST.ppt									# Conductivity file
InputFileSto=ANISOFLOW_TEST.ppt									# Conductivity file
InputFilePptByZones=ANISOFLOW_TEST.cvt 						# Conductivity by zones file
InputFileBC=ANISOFLOW_TEST.bc 									 	# Boundary Condition file
InputFileInitSol=ANISOFLOW_Sol.h5

OutputType=3 												# Output files type
OutputDir=out/ 												# Output directory

Scheme=1 													# Scheme type used 
Time=1														# Transitory boolean


# declare -a KSP=("cg" "gmres") #\

# declare -a KSP=("richardson" "chebyshev" "cg" "groppcg" "pipecg" 	\
	# "cgne" "nash" "stcg" "gltr" "fcg" "gmres" "fgmres" "lgmres"	\
	# "dgmres" "pgmres" "tcqmr" "bcgs" "ibcgs" "fbcgs" "fbcgsr"	\
	# "bcgsl" "cgs" "tfqmr" "cr" "pipecr" "lsqr" "preonly" "qcg"	\
	# "bicg" "minres" "symmlq" "lcd" "python" "gcr")

# declare -a PC=("none" "jacobi" "bjacobi")

# declare -a PC=("none" "jacobi" "sor" "lu" "shell" "bjacobi" 		\
# 	"mg" "eisenstat" "ilu" "icc" "asm" "gasm" "ksp" "composite" 	\
# 	"redundant" "spai" "nn" "cholesky" "pbjacobi" "mat" "hypre" 	\
# 	"parms" "fieldsplit" "tfs" "ml" "galerkin" "exotic" "cp" "bfbt"	\
# 	"lsc" "python" "pfmg" "syspfmg" "redistribute" "svd" "gamg"		\
# 	"sacusp" "sacusppoly" "bicgstabcusp" "ainvcusp" "bddc"			\ 
# 	"kaczmarz")

np=2 # for np in `seq 1 $np_max`; do
ksptype=gmres # 	for ksptype in "${KSP[@]}"; do
pctype=bjacobi # 		for pctype in "${PC[@]}"; do
			echo "*************************************************"
			echo "*************** Beginning new run ***************"
			echo "Options"
			echo "np: $np"
			echo "Solver: $ksptype"
			echo "Preconditioner: $pctype"
			echo "************** Initializing program *************"
			mpiexec -n $np ./$ProgramDir$ProgramName 			\
			-Input_type 				$InputType 				\
			-Input_type_tplgy 			$InputTypeTplgy 		\
			-Input_type_gmtry 			$InputTypeGmtry 		\
			-Input_type_bc	 			$InputTypeBC 	 		\
			-Input_type_init_sol	 	$InputTypeInitSol		\
			-Input_dir 					$InputDir 				\
			-Input_file_gmtry 			$InputFileGmtry 		\
			-Input_file_tplgy 			$InputFileTplgy 		\
			-Input_file_cvt 			$InputFileCvt 			\
			-Input_file_sto 			$InputFileSto 			\
			-Input_file_ppt_by_zones 	$InputFilePptByZones 	\
			-Input_file_bc 				$InputFileBC 			\
			-Input_file_init_sol		$InputFileInitSol		\
			-Output_type 				$OutputType 			\
			-Output_dir 				$OutputDir 				\
			-Run_options_scheme 		$Scheme 				\
			-Run_options_time 			$Time 					\
			-ksp_type 					$ksptype 				\
			-pc_type $pctype 				\
			-pc_sor_omega 1 -pc_sor_symmetric -pc_sor_its 1\
			-ksp_initial_guess_nonzero 							\
			-ksp_monitor 										\
			-ksp_monitor_lg_residualnorm 1 		\
			-log_view 					ascii:$OutputDir/ANISOFLOW-NP$np-KSP$ksptype-PC$pctype.log \
			| tee $OutputDir/ANISOFLOW-NP$np-KSP$ksptype-PC$pctype.log_term
			echo "***************** Ending program ****************"
# 		done
# 	done
# done
