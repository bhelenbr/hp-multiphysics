#!/bin/bash
#  Runs a thermal simulation of a NAND gate
#  

# This script is a fancy way of finding executables
# and is flexible enough to run on different platforms
# there are two main inputs: NP which is the number
# of processors to use (set to 0 for serial)
# and USE_VALGRIND which turns on debugging with valgrind (set to 1)
# if run from sun grid engine, NP will be set equal to NSLOTS which
# is the number of processors requested from the queue

INPT="run_b0.inpt"

# Set number of processors (overridden if using grid engine)
let NP=3

# Turn valgrind on/off
let USE_VALGRIND=0

# Uncomment to add petsc debug flags
#PETSC_FLAGS+=" -info -log_summary -ac-log_summary -memory_info -malloc_log -malloc_info -malloc_debug"
#PETSC_FLAGS+=" -fp_trap"
#PETSC_FLAGS+=" -on_error_attach_debugger gdb"
#PETSC_FLAGS+=" -start_in_debugger gdb"
#PETSC_FLAGS+=" -stop_for_debugger"

# Uncomment to set valgrind debugging parameters
VALGRIND_FLAGS+=" --track-origins=yes"
#VALGRIND_FLAGS+=" --leak-check=yes"
VALGRIND_FLAGS+=" --dsymutil=yes"
 
 
############################################
# This part should "just work" (in theory)
# Skip to end to modify what is actually done
############################################
# Detect if running from Sun Grid Engine
if [ -n "${NSLOTS}" ]; then
	# Running from gridengine
	let NP=${NSLOTS}
    MF="-machinefile ${TMPDIR}/machines"
fi

# For platform specific valgrind stuff
if [ ${NP} -gt 0 ]; then
	EXECNAME="tet_hp_mpi"
	if [ -e ${HOME}/Packages/bin/mpiexec ]; then
		export LD_PRELOAD="${HOME}/Packages/lib/valgrind/libmpiwrap-amd64-darwin.so"
		VALGRIND_FLAGS+=" --suppressions=${HOME}/Packages/share/openmpi/openmpi-valgrind.supp"
	else
		VALGRIND_FLAGS+=" --suppressions=/usr/local/src/mvapich2-1.9a/src/mpid/ch3/channels/mrail/src/hwloc/contrib/hwloc-valgrind.supp"
	fi
	EXECENV="mpiexec -np ${NP} ${MF}"
	INPTFILE="comm.inpt"
else
	EXECNAME="tet_hp"
	INPTFILE=${INPT}
fi


VALGRIND="valgrind  ${VALGRIND_FLAGS}"
MYPATH="${HOME}/bin"

#Set executable command
if [ ${USE_VALGRIND} -eq 0 ]; then
	HP="${EXECENV} ${MYPATH}/${EXECNAME} ${INPTFILE} ${PETSC_FLAGS}"
else
	HP=${EXECENV} ${VALGRIND} ${MYPATH}/${EXECNAME} ${INPTFILE} ${PETSC_FLAGS}
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

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm *

cp ../Inputs/* .

MESHNAME=$(mod_map -e run_b0.inpt b0_mesh)
tet_mesh -gl ${MESHNAME}

if [ ${NP} -gt 0 ]; then
	mod_map run_b0.inpt nblock ${NP}
	${HOME}/Codes/tet_mesh/partition2.bash ${INPT}
	sed 's/partition$/comm/g' partition.inpt > comm.inpt
fi


# Run the executable 
${HP}

cd ..
opendiff Baseline Results
