

ProgramDir=../../src/
ProgramName=ANISOFLOW

np=2														# Number of processors to be used


InputType=1													# Input files type
InputTypeTplgy=1 											# Input toplogy file type
InputTypeGmtry=2 											# Input gemoetry file type
InputTypeCvt=2
InputDir=in/ 												# Input directory
InputFileGmtry=sanpck.domnRST 											# Geometry file
InputFileTplgy=sanpck.acte 									# Topology file
InputFileCvt=sanpck.cvte												# Conductivity file
#InputFileCvtByZones= 										# Conductivity by zones file
InputFileBC=sanpck.bce 												# Boundary Condition file

OutputType=2 												# Output files type
OutputDir=out/ 												# Output directory

Scheme=1

mpiexec -n $np ./$ProgramDir$ProgramName -Input_type $InputType -Input_type_tplgy $InputTypeTplgy -Input_type_gmtry $InputTypeGmtry -Input_type_cvt $InputTypeCvt -Input_Dir $InputDir -Input_file_gmtry $InputFileGmtry -Input_file_cvt $InputFileCvt -Input_file_tplgy $InputFileTplgy -Input_file_bc $InputFileBC -Output_type $OutputType -Output_dir $OutputDir -Run_options_scheme $Scheme 
