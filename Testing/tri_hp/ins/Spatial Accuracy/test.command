#!/bin/bash
# Runs a Kovasznay solution problem on a series of meshes
# Checks spatial convergence rate
# Runs once with a square domain
# Runs second time with curved domain
# check Error.pdf & Curv_Error.pdf & rates.pdf

# This script is a fancy way of finding executables
# and is flexible enough to run on different platforms
# there are two main inputs: NP which is the number
# of processors to use (set to 0 for serial)
# and USE_VALGRIND which turns on debugging with valgrind (set to 1)
# if run from sun grid engine, NP will be set equal to NSLOTS which
# is the number of processors requested from the queue

# Set number of processors (overridden if using grid engine)
# For serial job set NP=0
let NP=1

# Turn valgrind on/off
let USE_VALGRIND=0

# Set executable to run
# Choose this for multigrid
#EXECNAME="tri_hp"
# Choose this for petsc
EXECNAME="tri_hp_petsc"

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


cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/*/*/*/*}/bin
export PATH=${PATH}:${BINDIR}

if [ -e Results ]; then
	cd Results
	rm *
else
	mkdir Results
	cd Results
fi

cp ../Inputs/* .

let nrefinements=8

let log2p=0
while [ $log2p -lt 3 ]; do
	let ngrid=2
	while [ $ngrid -lt $nrefinements ]; do
		mod_map run.inpt b0_mesh square${ngrid}.grd
		mod_map run.inpt log2p ${log2p}
#		mod_map run.inpt extra_finest_levels ${log2p}
#		mod_map run.inpt ngrid ${ngrid}
		mod_map run.inpt logfile output${ngrid}_p${log2p}
		${HP} run
		tail -2 output${ngrid}_p${log2p}_b0.log | head -1 | cut -d\  -f2,4,6,8,10,12 >> cnvg${log2p}.dat
		let ngrid=${ngrid}+1
	done
	let log2p=${log2p}+1
done



mod_map run.inpt A 0.2
cp run.inpt square2_bdry.inpt
tri_mesh -r square2 square3
tri_mesh -f square3 square3

let log2p=0
while [ $log2p -lt 3 ]; do
	let ngrid=3
	while [ $ngrid -lt $nrefinements ]; do
		mod_map run.inpt b0_mesh ./square${ngrid}.grd
		mod_map run.inpt log2p ${log2p}
#		mod_map run.inpt extra_finest_levels ${log2p}
#		mod_map run.inpt ngrid ${ngrid}
		mod_map run.inpt logfile curved${ngrid}_p${log2p}
		${HP} run
		tail -2 curved${ngrid}_p${log2p}_b0.log | head -1 | cut -d\  -f2,4,6,8,10,12 >> cnvg_curved${log2p}.dat
		
		let ngp=${ngrid}+1
		cp run.inpt square${ngrid}_bdry.inpt
		tri_mesh -r square${ngrid} square${ngp}
		tri_mesh -f square${ngp} square${ngp}
		let ngrid=${ngp}
	done
	let log2p=${log2p}+1
done

cd ..

./make_plot.command > Results/rates.dat

opendiff Results/ Baseline/
