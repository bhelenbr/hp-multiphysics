#!/bin/bash
#  This is a test of a mesh generation, coarsening, and translation
#  Works well with Laplacian deformation
#  Biharmonic deformation diverges?

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

# generate mesh and remove unnecessary data files
${EX} generate.inpt
mv rstrt1_b0.grd generate_b0.grd
NVRTX0=$(head -1 generate_b0.grd | cut -d" " -f 2)
COUNTER=0
while [  $COUNTER -lt $NVRTX0 ]; do
 echo "0.25" >> generate_b0.lngth
 let COUNTER=COUNTER+1 
done

mv rstrt1_b1.grd generate_b1.grd
NVRTX0=$(head -1 generate_b1.grd | cut -d" " -f 2)
COUNTER=0
while [  $COUNTER -lt $NVRTX0 ]; do
 echo "0.25" >> generate_b1.lngth
 let COUNTER=COUNTER+1 
done

${EX} coarsen.inpt

mv rstrt1_b0.grd coarsen_b0.grd
mv rstrt1_b1.grd coarsen_b1.grd
${EX} translate.inpt

# use opendiff (on OS X) to compare to a Baseline run
# can change to diff to do this on linux
cd ..
opendiff Baseline/ Results/
