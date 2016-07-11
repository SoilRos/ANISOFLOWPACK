

ProgramDir=../../src/
ProgramName=ANISOFLOW

np=1														# Number of processors to be used

InputType=1													# Input files type
InputTypeTplgy=1 											# Input toplogy file type
InputTypeGmtry=2 											# Input gemoetry file type
InputDir=in/ 												# Input directory
InputFileGmtry=Blesset_Grid 								# Geometry file
InputFileTplgy=Blesset_Tplgy 								# Topology file
InputFileCvt=Blesset_Cvt									# Conductivity file
InputFileCvtByZones=Blesset_CvtByZone 						# Conductivity by zones file
InputFileBC=Blesset_BC 									 	# Boundary Condition file

OutputType=3 												# Output files type
OutputDir=out/ 												# Output directory

Scheme=2 													# Scheme type used 
Time=0														# Transitory boolean

mpiexec -n $np ./$ProgramDir$ProgramName 			\
-Input_type 				$InputType 				\
-Input_type_tplgy 			$InputTypeTplgy 		\
-Input_type_gmtry 			$InputTypeGmtry 		\
-Input_dir 					$InputDir 				\
-Input_file_gmtry 			$InputFileGmtry 		\
-Input_file_tplgy 			$InputFileTplgy 		\
-Input_file_cvt 			$InputFileCvt 			\
-Input_file_cvt_by_zones 	$InputFileCvtByZones 	\
-Input_file_bc 				$InputFileBC 			\
-Output_type 				$OutputType 			\
-Output_dir 				$OutputDir 				\
-Run_options_scheme 		$Scheme 				\
-Run_options_time 			$Time 					\
-log_view 					ascii:ANISOFLOW_$np.log \
-ksp_monitor_lg_residualnorm 1