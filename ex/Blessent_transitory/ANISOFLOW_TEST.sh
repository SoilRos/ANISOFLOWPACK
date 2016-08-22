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

InputOutputVariables=										\
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
-Run_options_time 			$Time 


declare -a KSP=("richardson" "chebyshev" "cg" "groppcg" "pipecg") #\
	# "cgne" "nash" "stcg" "gltr" "fcg" "gmres" "fgmres" "lgmres"	\
	# "dgmres" "pgmres" "tcqmr" "bcgs" "ibcgs" "fbcgs" "fbcgsr"	\
	# "bcgsl" "cgs" "tfqmr" "cr" "pipecr" "lsqr" "preonly" "qcg"	\
	# "bicg" "minres" "symmlq" "lcd" "python" "gcr")

declare -a PC=("none" "jacobi" "sor" "lu" "shell" "bjacobi" 		\
	"mg" "eisenstat" "ilu" "icc" "asm" "gasm" "ksp" "composite" 	\
	"redundant" "spai" "nn" "cholesky" "pbjacobi" "mat" "hypre" 	\
	"parms" "fieldsplit" "tfs" "ml" "galerkin" "exotic" "cp" "bfbt"	\
	"lsc" "python" "pfmg" "syspfmg" "redistribute" "svd" "gamg"		\
	"sacusp" "sacusppoly" "bicgstabcusp" "ainvcusp" "bddc"			\ 
	"kaczmarz")

ksptype=cg
np=2														# Number of processors to be used

# for np in 1 2 4; do
	# for ksptype in "${KSP[@]}"; do
		for pctype in "${PC[@]}"; do
			if [ "$pctype" = "ilu" ]; then
				for level in 0 1 2; do
					echo "***** Beginning new run *****"
					echo mpiexec -n $np ./$ProgramDir$ProgramName 			\
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
					-pc_type 					$pctype 				\
					-pc_factor_levels 			$level					\
					-ksp_monitor 										\
					-log_view 					ascii:ANISOFLOW_NP$np-KSP$ksptype-PC$pctype-LVL$level.log #\
					# >> ANISOFLOW_NP$np_KSP$ksptype_PC$pctype_LVL$level.log_term
				done
			else
				echo "***** Beginning new run *****"
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
				-pc_type 					$pctype 				\
				-ksp_monitor 										\
				-log_view 					ascii:ANISOFLOW_NP$np-KSP$ksptype-PC$pctype.log #\
				# >> ANISOFLOW_NP$np_KSP$ksptype_PC$pctype_LVL$level.log_term
			fi
		done
	# done
# done



# mpiexec -n $np ./$ProgramDir$ProgramName 			\
# -Input_type 				$InputType 				\
# -Input_type_tplgy 			$InputTypeTplgy 		\
# -Input_type_gmtry 			$InputTypeGmtry 		\
# -Input_type_bc	 			$InputTypeBC 	 		\
# -Input_type_init_sol	 	$InputTypeInitSol		\
# -Input_dir 					$InputDir 				\
# -Input_file_gmtry 			$InputFileGmtry 		\
# -Input_file_tplgy 			$InputFileTplgy 		\
# -Input_file_cvt 			$InputFileCvt 			\
# -Input_file_sto 			$InputFileSto 			\
# -Input_file_ppt_by_zones 	$InputFilePptByZones 	\
# -Input_file_bc 				$InputFileBC 			\
# -Input_file_init_sol		$InputFileInitSol		\
# -Output_type 				$OutputType 			\
# -Output_dir 				$OutputDir 				\
# -Run_options_scheme 		$Scheme 				\
# -Run_options_time 			$Time 					\
# -log_view 					ascii:ANISOFLOW_$np.log 
#-ksp_monitor_lg_residualnorm 1
