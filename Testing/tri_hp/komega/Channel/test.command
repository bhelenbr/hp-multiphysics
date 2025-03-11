#!/bin/bash

# Gird convergence study of fully-developed channel flow driven by gravity using periodic boundary conditions. The grids
# correspond to ones studied with the 1-D MATLAB code. The solution from the 1-D MATLAB code is used to restart 
# tri_hp_petsc. The solutions for both Menter's and Wilcox's BCs are obtained with p = 1, 2, and 4

# Need to use DIRK1, and set susk = susomg=0. 
# Can choose between CALC_TAU1, CALC_TAU2 and WILCOX_1988, WILCOX_1988KL and WILCOX2006.


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
HP="mpiexec -np ${NP} tri_hp_petsc run.inpt ${PETSC_FLAGS}"

# Make Results directory
if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm *

# copy input files into results directory
cp ../Inputs/* .

# generate mesh 
tri_mesh generate

# Run the executable (for each value of log2p, the Wilcox's BC is done then Menter's)

${HP}


if [ "$?" -eq "0" ]; then
	let RESTART=$(grep TIMESTEP: out_b0.log | tail -1 | cut -d\  -f2)
	let NC=9
	mod_map run.inpt nconv ${NC}
	mod_map run.inpt restart ${RESTART}
	${HP}

	while [ ${NC} -lt 32 ]; do		

		let RESTART=${RESTART}+10
		let NC=${NC}+1
		mod_map run.inpt nconv ${NC}
		mod_map run.inpt restart ${RESTART}
		${HP}

		while [ "$?" -ne "0" ]; do
	 		rm core*
	 		rm neg*
	 		rm abort*
			let NC=${NC}-1
			let RESTART=${RESTART}-1
			mod_map run.inpt restart ${RESTART}
			mod_map run.inpt nconv ${NC}
			${HP}
		done
	done
	let RESTART=${RESTART}+10
	mod_map run.inpt dtinv 0
	mod_map run.inpt restart ${RESTART}
	${HP}
	es=$?
fi

if [ "${es}" -eq "0" ]; then
	mod_map run.inpt log2p 1
	let RESTART=$(grep TIMESTEP: out_b0.log | tail -1 | cut -d\  -f2)
	mod_map run.inpt restart ${RESTART}
	${HP}
	es=$?
fi

if [ "${es}" -eq "0" ]; then
	mod_map run.inpt log2p 2
	let RESTART=$(grep TIMESTEP: out_b0.log | tail -1 | cut -d\  -f2)
	mod_map run.inpt restart ${RESTART}
	${HP}
fi


