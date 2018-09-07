#!/bin/bash
#  This is a test of a mesh generation, coarsening, and translation
#  Works well with Laplacian deformation
#  Biharmonic deformation diverges?

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
