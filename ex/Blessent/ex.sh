#!/bin/bash

max_np=12
#PBS -N ANISOFLOW			# nombre identificador del trabajo
#PBS -q cola$max_np			  	 # Nombre de la cola de ejecucion
#PBS -l nodes=1:ppn=$max_np           # nodes=1 numero de nodos a emplear
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
InputFileGmtry=Blesset_Grid 								# Geometry file
InputFileTplgy=Blesset_Tplgy 								# Topology file
InputFileCvt=Blesset_Cvt									# Conductivity file
InputFileSto=Blesset_Cvt									# Conductivity file
InputFilePptByZones=Blesset_CvtByZone 						# Conductivity by zones file
InputFileBC=Blesset_BC 									 	# Boundary Condition file
InputFileInitSol=ANISOFLOW_Sol.h5

OutputType=3 												# Output files type
OutputDir=out/ 												# Output directory

Scheme=1 													# Scheme type used 
Time=0														# Transitory boolean

np=$max_np

mpiexec -n $np ./$ProgramDir$ProgramName 			\
-Input_type 				$InputType 				\
-Input_type_tplgy 			$InputTypeTplgy 		\
-Input_type_gmtry 			$InputTypeGmtry 		\
-Input_type_bc	 			$InputTypeBC 	 		\
-Input_dir 					$InputDir 				\
-Input_file_gmtry 			$InputFileGmtry 		\
-Input_file_tplgy 			$InputFileTplgy 		\
-Input_file_cvt 			$InputFileCvt 			\
-Input_file_sto 			$InputFileSto 			\
-Input_file_ppt_by_zones 	$InputFilePptByZones 	\
-Input_file_bc 				$InputFileBC 			\
-Output_type 				$OutputType 			\
-Output_dir 				$OutputDir 				\
-Run_options_scheme 		$Scheme 				\
-Run_options_time 			$Time 					\
-log_view 					ascii:ANISOFLOW_$np.log \
-ksp_monitor_lg_residualnorm 1  					\
-ksp_initial_guess_nonzero -ksp_monitor -ksp_type cg -pc_type pbjacobi
