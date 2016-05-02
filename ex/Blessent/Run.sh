

ProgramDir=../../src/
ProgramName=ANISOFLOW.exe

np=1														# Number of processors to be used


InputType=1													# Input files type
InputDir=in/ 												# Input directory
InputFileGmtry=tsim_USMH.asc 								# Geometry file
InputFileCvt=matrix.mprops									# Conductivity file
InputFileCvtByZones=tsim_USMH.asc 							# Conductivity by zones file
InputFileBC=grid_400_400.nch_nprop_list.lateral_boundary 	# Boundary Condition file

Scheme=2

mpiexec -n $np ./$ProgramDir$ProgramName -Input_type $InputType -Input_Dir $InputDir -Input_file_gmtry $InputFileGmtry -Input_file_cvt $InputFileCvt -Input_file_cvt_by_zones $InputFileCvtByZones -Input_file_bc $InputFileBC -Run_options_scheme $Scheme
