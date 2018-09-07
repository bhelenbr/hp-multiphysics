#!/bin/bash
#  Runs a deforming adapting mesh for a swimming flagellum
#  Tests the mesh deformation & adaption routines
#  Also tests mesh coarsening.  
#  When using OLDRECONNECT triangulate fails because boundary sides cross on coarse meshes
#  Works with coarsen2.

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
EXECNAME="tri_mesh"

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

${EX} generate

cp generate.inpt translate.inpt
mod_map translate.inpt b0_mesh rstrt1_b0.grd
mod_map translate.inpt ngrid 3
mod_map translate.inpt ntstep 20
mod_map translate.inpt dtinv 20.0
mod_map translate.inpt ncycle 100
mod_map translate.inpt logfile translate
mod_map translate.inpt vwcycle 2
${EX} translate.inpt

# use opendiff (on OS X) to compare to a Baseline run
# can change to diff to do this on linux
cd ..
opendiff Baseline/ Results/