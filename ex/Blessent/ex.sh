

ProgramDir=../../src/
ProgramName=ANISOFLOW

np=1														# Number of processors to be used

InputType=1													# Input files type
InputTypeTplgy=1 											# Input toplogy file type
InputTypeGmtry=2 											# Input gemoetry file type
InputTypeBC=1 												# Input boundary conditions file type
InputTypeInitSol=3
InputDir=in/ 												# Input directory
InputFileGmtry=Blesset_Grid 								# Geometry file
InputFileTplgy=Blesset_Tplgy 								# Topology file
InputFileCvt=Blesset_Cvt									# Conductivity file
InputFileSto=Blesset_Cvt									# Conductivity file
InputFilePptByZones=Blesset_CvtByZone 						# Conductivity by zones file
InputFileBC=Blesset_BC 									 	# Boundary Condition file
InputFileInitSol=ANISOFLOW_Sol.h5

OutputType=3 												# Output files type
OutputDir=out/ 												# Output directory

Scheme=2 													# Scheme type used 
Time=1														# Transitory boolean

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
-log_view 					ascii:ANISOFLOW_$np.log \
-ksp_monitor_lg_residualnorm 1