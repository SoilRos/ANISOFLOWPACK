 # Execution Line must be: ./25x25x25_PerformanceTest.sh np_max
 CantTest=$1
 np=1
 declare -a KSP=("cg" "gmres")
 declare -a PC=("none" "jacobi" "bjacobi" "mg")
 for ContTest in `seq 1 ${CantTest}`; do
     np=`expr 2 \* ${np}
     for ksptype in "${KSP[@]}"; do
         for pctype in "${PC[@]}"; do
             echo "*************************************************"
             echo "*************** Beginning new run ***************"
             echo "Options"
             echo "np: $np"
             echo "Solver: $ksptype"
             echo "Preconditioner: $pctype"
             echo "************** Initializing program *************"
             mpiexec -n ${np} \
             ../../src/ANISOFLOW \
             -Input_type_gmtry 2 \
             -Input_file_gmtry 25x25x25_gmtry.txt \
             -Input_type_tplgy 0 \
             -Input_file_ppt 25x25x25_properties.txt \
             -Input_file_ppt_by_zones 25x25x25_zones.txt \
             -Input_file_bc 25x25x25_bc.txt \
             -Input_type_bc 2 \
             -Run_options_scheme 2 \
             -Output_type_ppt 3 \
             -Project_name 25x25x25\
             -ksp_type ${ksptype} \
             -pc_type  ${pctype} \
              
             mpiexec -n ${np} \
             ../../src/ANISOFLOW \
             -Input_type_gmtry 2 \
             -Input_file_gmtry 25x25x25_gmtry.txt \
             -Input_type_tplgy 0 \
             -Input_file_ppt 25x25x25_properties.txt \
             -Input_file_ppt_by_zones 25x25x25_zones.txt \
             -Input_file_bc 25x25x25_bc.txt \
             -Input_type_bc 2 \
             -Run_options_scheme 2 \
             -Run_options_time 1 \
             -Output_type_ppt 3 \
             -Project_name 25x25x25 \
             -Input_type_init_sol 3 \
             -Input_file_init_sol 25x25x25_Sol.h5 \
             -ksp_type ${ksptype} \
             -pc_type  ${pctype} \
             -log_view ascii:25x25x25_NP${np}_KSP${ksptype}_PC${pctype}.log \
             | tee 25x25x25_NP${np}_KSP${ksptype}_PC${pctype}.log_term 
             echo "***************** Ending program ****************"
         done
     done
 done
