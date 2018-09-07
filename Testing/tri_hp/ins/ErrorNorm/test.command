#!/bin/bash
#  Checks error measure effect on convergence and mesh adaptation
#  Uses Kovaznay solution and decrease error tolerance
#  Then calculates L2 error as a function of error tolerance

# This script is a fancy way of finding executables
# and is flexible enough to run on different platforms
# there are two main inputs: NP which is the number
# of processors to use (set to 0 for serial)
# and USE_VALGRIND which turns on debugging with valgrind (set to 1)
# if run from sun grid engine, NP will be set equal to NSLOTS which
# is the number of processors requested from the queue

# Set number of processors (overridden if using grid engine)
# For serial job set NP=0
let NP=0

# Turn valgrind on/off
let USE_VALGRIND=0

# Set executable to run
# Choose this for multigrid
#EXECNAME="tri_hp"
# Choose this for petsc
EXECNAME="tri_hp"

# Uncomment to add petsc debug flags
#PETSC_FLAGS+=" -info -log_summary -ac-log_summary -memory_info -malloc_log -malloc_info -malloc_debug"
#PETSC_FLAGS+=" -fp_trap"
#PETSC_FLAGS+=" -on_error_attach_debugger gdb"
#PETSC_FLAGS+=" -start_in_debugger gdb"
#PETSC_FLAGS+=" -stop_for_debugger"

# Uncomment to set valgrind debugging parameters
VALGRIND_FLAGS+=" --track-origins=yes"
VALGRIND_FLAGS+=" --leak-check=full"
VALGRIND_FLAGS+=" --dsymutil=yes"
 
 
############################################
# This part should "just work" (in theory)
# Skip to end to modify what is actually done
############################################
# Detect if running from Sun Grid Engine
if [ "${PE}" = "mpich" ]; then
	# Running from gridengine
	let NP=${NSLOTS}
    MF="-machinefile ${TMPDIR}/machines"
fi

# For platform specific valgrind stuff
if [ "${NP}" -gt 0 ]; then
	if [ -e "${HOME}/Packages/lib/valgrind/libmpiwrap-amd64-darwin.so" ]; then
		export LD_PRELOAD="${HOME}/Packages/lib/valgrind/libmpiwrap-amd64-darwin.so"
		VALGRIND_FLAGS+=" --suppressions=${HOME}/Packages/share/openmpi/openmpi-valgrind.supp"
	else
		VALGRIND_FLAGS+=" --suppressions=/usr/local/src/mvapich2-1.9a/src/mpid/ch3/channels/mrail/src/hwloc/contrib/hwloc-valgrind.supp"
	fi
	EXECENV="mpiexec -np ${NP} ${MF}"
fi

VALGRIND="valgrind  ${VALGRIND_FLAGS}"

#Set executable command
if [ "${USE_VALGRIND}" -eq 0 ]; then
	HP="${EXECENV} ${EXECNAME} run.inpt ${PETSC_FLAGS}"
else
	HP=${EXECENV} ${VALGRIND} ${MYPATH}/${EXECNAME} run.inpt ${PETSC_FLAGS}
fi

#############################################
# End of automated part HP is now the correct
# string needed to run the executable     
# can just copy and paste above part into shell if you want
# to run the executable without using this script
# Paste above into shell to define HP
# then type ${HP} to run
#############################################


# cd to the directory where the script resides.
# Necessary for platforms where script can be 
# double clicked from gui
cd "$(dirname "$0")"

# Make Results directory
if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm -rf *

# copy input files into results directory
cp ../Inputs/* .

let log2p=1
while [ $log2p -lt 3 ]; do
	let nerror=0
	while [ $nerror -lt 3 ]; do
		NAME=output${nerror}_p${log2p}
		mkdir ${NAME}
		cd ${NAME}
		cp ../run.inpt .
		mod_map run.inpt b0_mesh ${HOME}/Codes/Testing/grids/SQUARE/INOUT/square2
		mod_map run.inpt log2p ${log2p}
		mod_map run.inpt extra_finest_levels ${log2p}
		mod_map run.inpt logfile ${NAME}
		mod_map run.inpt nerror ${nerror}
		let ngrid=${nerror}+1
		mod_map run.inpt ngrid ${ngrid}

		${HP} run
		grep error_target ${NAME}_b0.log | tail -1 | cut -d\  -f3 | tr -d '\n' >> ../cnvg${log2p}.dat
		echo -n ' ' >> ../cnvg${log2p}.dat
		grep DOF ${NAME}_b0.log | tail -1 | cut -d\  -f6 | tr -d '\n' >> ../cnvg${log2p}.dat
		echo -n ' ' >> ../cnvg${log2p}.dat
		grep L_2 ${NAME}_b0.log | tail -1 | cut -d\  -f2,4,6,8,10,12 >> ../cnvg${log2p}.dat
		let nerror=${nerror}+1
		cd ..
	done
	let log2p=${log2p}+1
done

mod_map run.inpt A 0.05
let log2p=1
while [ $log2p -lt 3 ]; do
	let nerror=0
	while [ $nerror -lt 3 ]; do
		NAME=curved${nerror}_p${log2p}
		mkdir ${NAME}
		cd ${NAME}
		cp ../run.inpt .
		mod_map run.inpt b0_mesh ${HOME}/Codes/Testing/grids/SQUARE/INOUT/square2
		mod_map run.inpt log2p ${log2p}
		mod_map run.inpt extra_finest_levels ${log2p}
		mod_map run.inpt logfile ${NAME}
		mod_map run.inpt nerror ${nerror}
		let ngrid=${nerror}+1
		mod_map run.inpt ngrid ${ngrid}
		${HP} run
		grep error_target ${NAME}_b0.log | tail -1 | cut -d\  -f3 | tr -d '\n' >> ../cnvg_curved${log2p}.dat
		echo -n ' ' >> ../cnvg_curved${log2p}.dat
		grep DOF ${NAME}_b0.log | tail -1 | cut -d\  -f6 | tr -d '\n' >> ../cnvg_curved${log2p}.dat
		echo -n ' ' >> ../cnvg_curved${log2p}.dat
		grep L_2 ${NAME}_b0.log | tail -1 | cut -d\  -f2,4,6,8,10,12 >> ../cnvg_curved${log2p}.dat
		let nerror=${nerror}+1
		cd ..
	done
	let log2p=${log2p}+1
done

cd ..

opendiff Results/ BASELINE/
