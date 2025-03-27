#!/bin/bash

# Select the version 

version=1988
# version=1988KL
# version=2006
# version=2006NSL



# Number of processors
let NP=1

# cd to the directory where the script resides.
cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/*/*/*/*}/bin
export PATH=${PATH}:${BINDIR}

# Uncomment to add petsc debug flags
#PETSC_FLAGS+=" -info -log_summary -ac-log_summary -memory_info -malloc_log -malloc_info -malloc_debug"
#PETSC_FLAGS+=" -fp_trap"
#PETSC_FLAGS+=" -on_error_attach_debugger gdb"
#PETSC_FLAGS+=" -start_in_debugger gdb"
#PETSC_FLAGS+=" -stop_for_debugger"

# Set executable command
HP="mpiexec -np ${NP} tri_hp_petsc" 
HPAXI="mpiexec -np ${NP} tri_hp_axi_petsc" 


# Make Results directory
if [ -e Results2D${version} ]; then
	cd Results2D${version}
else
	mkdir Results2D${version}
	cd Results2D${version}
fi
rm *

# copy input files into results directory
cp ../Inputs/* .

# generate mesh 
tri_mesh generate


${HP} run2D${version}.inpt ${PETSC_FLAGS}

mod_map -c run2D${version}.inpt rsdl_debug
mv rstrt1_b0.grd square1_b0.grd

let nrefinements=3
let log2p=2
while [ $log2p -lt 3 ]; do
	let ngrid=1
	while [ $ngrid -le $nrefinements ]; do
		mod_map run2D${version}.inpt b0_mesh square${ngrid}_b0.grd
		mod_map run2D${version}.inpt log2p ${log2p}
		mod_map run2D${version}.inpt logfile output${ngrid}_p${log2p}
		${HP} run2D${version}.inpt ${PETSC}
		tail -2 output${ngrid}_p${log2p}_b0.log | head -1 | cut -d\  -f2,4,6,8,10,12,14,16,18,20 >> cnvg${log2p}.dat
		let ngp=${ngrid}+1
		tri_mesh -r square${ngrid}_b0.grd square${ngp}_b0.grd
		let ngrid=${ngp}
	done
	let log2p=${log2p}+1
done

cd ..

./make_plot.command ${version}

# ./make_plot.command > Results/rates.dat





# ${HP} run2D${version}.inpt ${PETSC_FLAGS}
# mod_map run2D${version}.inpt auto_timestep_tries 0
# mod_map run2D${version}.inpt dtinv 0.0
# mod_map run2D${version}.inpt ntstep 1
# mod_map run2D${version}.inpt restart 100
# mod_map run2D${version}.inpt ncycle 20
# # # mod_map run2D${version}.inpt petsc "-ksp_type gmres -ksp_max_it 1000 -ksp_atol 1e-10 -pc_type asm -pc_asm_subiterations 5"
# ${HP} run2D${version}.inpt ${PETSC_FLAGS}


# cd ..

# # Now run the axisymmetric case

# # Make Results directory
# if [ -e ResultsAxi${version} ]; then
# 	cd ResultsAxi${version}
# else
# 	mkdir ResultsAxi${version}
# 	cd ResultsAxi${version}
# fi
# rm *

# # copy input files into results directory
# cp ../Inputs/* .

# # generate mesh 
# tri_mesh generate

# ${HPAXI} runAxi${version}.inpt ${PETSC_FLAGS}



