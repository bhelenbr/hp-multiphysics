#!/bin/bash
#  Runs a deforming adapting mesh for a swimming flagellum
#  Tests the mesh deformation & adaption routines
#  Also tests mesh coarsening.  
#  When using OLDRECONNECT triangulate fails because boundary sides cross on coarse meshes
#  Works with coarsen2.

# cd to the directory where the script resides.
# Necessary for platforms where script can be 
# double clicked from gui
cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/*/*/*}/bin
echo ${BINDIR}
export PATH=${PATH}:${BINDIR}
EX="tri_mesh"

# Uncomment to set valgrind debugging parameters
VALGRIND_FLAGS+=" --track-origins=yes"
VALGRIND_FLAGS+=" --leak-check=full"
VALGRIND_FLAGS+=" --dsymutil=yes"
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