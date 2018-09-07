#!/bin/bash
# Test of partitioning script:
# partitions mesh then runs a deforming mesh iteration
# in parallel on the partitioned mesh

# This script is a fancy way of finding executables
# and is flexible enough to run on different platforms
# there are two main inputs: NP which is the number
# of processors to use (set to 0 for serial)
# and USE_VALGRIND which turns on debugging with valgrind (set to 1)
# if run from sun grid engine, NP will be set equal to NSLOTS which
# is the number of processors requested from the queue

# Turn valgrind on/off
let USE_VALGRIND=0

# Set executable to run
EXECNAME="tri_mesh_mpi"

# Uncomment to set valgrind debugging parameters
VALGRIND_FLAGS+=" --track-origins=yes"
VALGRIND_FLAGS+=" --leak-check=full"
VALGRIND_FLAGS+=" --dsymutil=yes"
VALGRIND="valgrind  ${VALGRIND_FLAGS}"

#Set executable command
if [ "${USE_VALGRIND}" -eq 0 ]; then
	EX="${EXECNAME}"
else
	EX=${VALGRIND} ${EXECNAME}
fi

# Detect if running from Sun Grid Engine
if [ "${PE}" = "mpich" ]; then
	# Running from gridengine
	let NP=${NSLOTS}
    MF="-machinefile ${TMPDIR}/machines"
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

# generate mesh and remove unnecessary data files
tri_mesh generate
mv rstrt1_b0.grd square_hole.grd

# BASELINE CASE
mkdir partition1
cd partition1
cp ../translate.inpt .
tri_mesh translate.inpt
cd ..

NPROC=(2 4 8)

let npc=0
while [ $npc -lt ${#NPROC[@]} ]; do
	let np=${NPROC[npc]}
	mkdir partition${np}
	cd partition${np}
	cp ../translate.inpt .
	mod_map translate.inpt nblock ${np}
	partition.bash translate.inpt
	mpiexec -np ${np} ${MF} ${EX} partition.inpt
	cd ..
	let npc=$npc+1
done

cd ..

opendiff Baseline/ Results/
