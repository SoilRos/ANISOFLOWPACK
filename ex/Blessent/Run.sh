

ProgramDir=../../src/
ProgramName=ANISOFLOW.exe

np=1														# Number of processors to be used


InputType=1													# Input files type
InputDir=in/ 												# Input directory
InputFileGmtry=tsim_USMH.asc 								# Geometry file
InputFileCvt=matrix.mprops									# Conductivity file
InputFileCvtByZones=tsim_USMH.asc 							# Conductivity by zones file
InputFileBC=BCtempDT 									 	# Boundary Condition file

OuputType=1 												# Ouput files type
OuputDir=out/ 												# Ouput directory

Scheme=2 													# Scheme type used 
Time=0														#

mpiexec -n $np ./$ProgramDir$ProgramName -Input_type $InputType -Input_dir $InputDir -Input_file_gmtry $InputFileGmtry -Input_file_cvt $InputFileCvt -Input_file_cvt_by_zones $InputFileCvtByZones -Input_file_bc $InputFileBC -Ouput_type $OuputType -Ouput_dir $OuputDir -Run_options_scheme $Scheme -Run_options_time $Time -log_view ascii:ANISOFLOW.log
