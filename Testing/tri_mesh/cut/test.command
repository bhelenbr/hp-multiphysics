#!/bin/bash
#  This is a test of the cut routine which cuts a mesh along an analytic curve
#  then smooth the boundary 

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
	EX="${EXECNAME} generate.inpt"
else
	EX=${VALGRIND} ${EXECNAME} generate.inpt
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

tri_mesh -z square32

# MAKE A NEW LENGTH INPUT FILE SO THAT MESH GETS COARSENED
NVRTX0=$(head -1 cut0.grd | cut -d" " -f 2)
COUNTER=0
while [  $COUNTER -lt $NVRTX0 ]; do
 echo "0.031" >> cut0.vlngth
 let COUNTER=COUNTER+1 
done

# MAKE A NEW LENGTH INPUT FILE SO THAT MESH GETS COARSENED
NVRTX1=$(head -1 cut1.grd | cut -d" " -f 2)
COUNTER=0
while [  $COUNTER -lt $NVRTX1 ]; do
 echo "0.031" >> cut1.vlngth
 let COUNTER=COUNTER+1 
done

${EX} generate

# use opendiff (on OS X) to compare to a Baseline run
# can change to diff to do this on linux
cd ..
opendiff Baseline/ Results/
