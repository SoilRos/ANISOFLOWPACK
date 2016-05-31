

ProgramDir=../../src/
ProgramName=ANISOFLOW.exe

np=1														# Number of processors to be used


InputType=1													# Input files type
InputTypeTplgy=2 											# Input toplogy file type
InputTypeGmtry=2 											# Input gemoetry file type
InputDir=in/ 												# Input directory
InputFileGmtry=Blesset_Grid 								# Geometry file
InputFileTplgy=Blesset_Tplgy 								# Topology file
InputFileCvt=Blesset_Cvt									# Conductivity file
InputFileCvtByZones=Blesset_CvtByZone 						# Conductivity by zones file
InputFileBC=Blesset_BC 									 	# Boundary Condition file

OuputType=2 												# Ouput files type
OuputDir=out/ 												# Ouput directory

Scheme=1 													# Scheme type used 
Time=0														#

mpiexec -n $np ./$ProgramDir$ProgramName -Input_type $InputType -Input_type_tplgy $InputTypeTplgy -Input_type_gmtry $InputTypeGmtry -Input_dir $InputDir -Input_file_gmtry $InputFileGmtry -Input_file_tplgy $InputFileTplgy -Input_file_cvt $InputFileCvt -Input_file_cvt_by_zones $InputFileCvtByZones -Input_file_bc $InputFileBC -Ouput_type $OuputType -Ouput_dir $OuputDir -Run_options_scheme $Scheme -Run_options_time $Time -log_view ascii:ANISOFLOW.log #-ksp_monitor_lg_residualnorm 1
