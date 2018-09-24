#!/bin/bash
# Runs a case with a solidifying surface 
# And deforms the solid surface so that it moves with the solid pull velocity
# not sure if this has worked yet.  Not sure why this test is tied in with
# the melting surface, but oh well.  The triple junction for this case is constrained to 
# move horizontally

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
#VALGRIND_FLAGS+=" --leak-check=yes"
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
if [ ${NP} -gt 0 ]; then
	if [ -e "${HOME}/Packages/lib/valgrind/libmpiwrap-amd64-darwin.so" ]; then
		export LD_PRELOAD="${HOME}/Packages/lib/valgrind/libmpiwrap-amd64-darwin.so"
		VALGRIND_FLAGS+=" --suppressions=${HOME}/Packages/share/openmpi/openmpi-valgrind.supp"
	else
		VALGRIND_FLAGS+=" --suppressions=/usr/local/src/mvapich2-1.9a/src/mpid/ch3/channels/mrail/src/hwloc/contrib/hwloc-valgrind.supp"
	fi
	EXECENV="mpiexec -np ${NP} ${MF}"
fi

VALGRIND="valgrind  ${VALGRIND_FLAGS}"
MYPATH="${HOME}/bin"

#Set executable command
if [ ${USE_VALGRIND} -eq 0 ]; then
	HP="${EXECENV} ${MYPATH}/${EXECNAME} run.inpt ${PETSC_FLAGS}"
else
	HP="${EXECENV} ${VALGRIND} ${MYPATH}/${EXECNAME} run.inpt ${PETSC_FLAGS}"
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

# These are the steps to run
# Get a reasonable steady temperature field
# Before starting the melting 
# mod_map -d  means delete that line from the input file
# mod_map run.inpt keyword value ## Changes the value in the input file
cp run.inpt startup.inpt
mod_map run.inpt b0_v1_type plain
mod_map run.inpt b0_v2_type plain
mod_map run.inpt b0_v3_type plain
mod_map run.inpt b0_v1_hp_type plain
mod_map run.inpt b0_v2_hp_type plain
mod_map run.inpt b0_v3_hp_type plain
mod_map run.inpt b1_v1_type plain
mod_map run.inpt b1_v2_type plain
mod_map run.inpt b1_v4_type plain
mod_map run.inpt b1_v1_hp_type plain
mod_map run.inpt b1_v2_hp_type plain
mod_map run.inpt b1_v4_hp_type plain
mod_map run.inpt b0_s3_hp_type inflow
mod_map run.inpt b0_s3_type symbolic
mod_map run.inpt b1_s3_type symbolic
mod_map run.inpt b1_s3_hp_type dirichlet
mod_map run.inpt b1_s7_type symbolic
mod_map run.inpt b1_s7_hp_type plain
mod_map run.inpt b1_s7_flux0 convflux+radiation*epss
mod_map run.inpt adapt 0
mod_map -c run.inpt mesh_movement
mod_map run.inpt ntstep 1
mod_map run.inpt dtinv 0.0
mod_map run.inpt restart_interval 1
${HP}
if [ "$?" -ne "0" ]; then
  echo "exited from initial T solution"
  exit 1
fi 
let RESTART=0
rename rstrt1 rstrt0 *.nc

mv startup.inpt run.inpt
# 100 steps with small time step
mod_map run.inpt ntstep 10
#mod_map run.inpt restart ${RESTART}
mod_map run.inpt restart_interval 100
${HP}
if [ "$?" -ne "0" ]; then
  echo "exited from unsteady evolution 1"
  exit 1
fi 

rm core*
rm abort*
rm neg*

cd ..

