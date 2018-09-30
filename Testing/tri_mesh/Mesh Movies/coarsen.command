#!/bin/bash

# This tests refinement and coarsening
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

if [ -e Results_Coarsen ]; then
	cd Results_Coarsen
else
	mkdir Results_Coarsen
	cd Results_Coarsen
fi
rm *

cp ../Inputs/* .

# MAKE A NEW LENGTH INPUT FILE SO THAT tri_mesh GETS COARSENED
NVRTX0=$(head -1 generate_b0.grd | cut -d" " -f 2)
COUNTER=0
while [  $COUNTER -lt $NVRTX0 ]; do
 echo "0.25" >> generate_b0.lngth
 let COUNTER=COUNTER+1 
done

echo "Beginning Coarsen"
${EX} coarsen
echo "Ending Coarsen"

cd ..

opendiff Results_Coarsen/ Results_Coarsen_Baseline/
