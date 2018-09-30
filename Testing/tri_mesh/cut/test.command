#!/bin/bash
#  This is a test of the cut routine which cuts a mesh along an analytic curve
#  then smooth the boundary 

# cd to the directory where the script resides.
# Necessary for platforms where script can be 
# double clicked from gui
cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/*/*/*}/bin
echo ${BINDIR}
export PATH=${PATH}:${BINDIR}

# Uncomment to set valgrind debugging parameters
VALGRIND_FLAGS+=" --track-origins=yes"
VALGRIND_FLAGS+=" --leak-check=full"
VALGRIND_FLAGS+=" --dsymutil=yes"

EX="tri_mesh"
#EX="valgrind  ${VALGRIND_FLAGS} tri_mesh"

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
