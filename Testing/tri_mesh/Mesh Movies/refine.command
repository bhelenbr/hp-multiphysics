#!/bin/bash

# This test refinement and coarsening
# with a non-uniform resolution
# Uncomment the DEBUG_ADAPT line
# in rebay.pdf to see the process for refinement
# Uncomment the DEBUG_ADAPT line
# in yaber.pdf to see the process for coarsening

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

if [ -e Results_Refine ]; then
	cd Results_Refine
else
	mkdir Results_Refine
	cd Results_Refine
fi
rm *

cp ../Inputs/* .
echo "Beginning Generate"
${EX} generate
echo "Ending Generate"
cd ..

opendiff Results_Refine/ Results_Refine_Baseline/
