

ProgramDir=../../src/
ProgramName=ANISOFLOW.exe

np=2														# Number of processors to be used


InputType=1													# Input files type
InputTypeTplgy=2 											# Input toplogy file type
InputTypeGmtry=2 											# Input gemoetry file type
InputDir=in/ 												# Input directory
InputFileGmtry=sanpck.domnRST 											# Geometry file
InputFileCvt=												# Conductivity file
InputFileCvtByZones= 										# Conductivity by zones file
InputFileBC= 												# Boundary Condition file

Scheme=1

mpiexec -n $np ./$ProgramDir$ProgramName -Input_type $InputType -Input_type_tplgy $InputTypeTplgy -Input_type_gmtry $InputTypeGmtry -Input_Dir $InputDir -Input_file_gmtry $InputFileGmtry -Input_file_cvt $InputFileCvt -Input_file_cvt_by_zones $InputFileCvtByZones -Input_file_bc $InputFileBC -Run_options_scheme $Scheme
