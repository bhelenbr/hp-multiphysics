#!/bin/bash
#  Put description of case here
#  

# This script is a fancy way of finding executables
# and is flexible enough to run on different platforms
# there are two main inputs: NP which is the number
# of processors to use (set to 0 for serial)
# and USE_VALGRIND which turns on debugging with valgrind (set to 1)
# if run from sun grid engine, NP will be set equal to NSLOTS which
# is the number of processors requested from the queue

# Set number of processors (overridden if using grid engine)
# For serial job set NP=0
let NP=2

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
rm *

# copy input files into results directory
cp ../Inputs/* .

# generate mesh and remove unnecessary data files
tri_mesh generate
rm data*.grd

# Run the executable 
${HP}

cd ..

# use opendiff (on OS X) to compare to a Baseline run
# can change to diff to do this on linux
opendiff Baseline/ Results/
