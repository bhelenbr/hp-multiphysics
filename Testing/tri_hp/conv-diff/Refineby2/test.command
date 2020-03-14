#!/bin/bash

# Tests adaptation of a function that should be represented exactly on the mesh = u = x^4
# coarsens the mesh and refines the mesh to make sure that is maintained.
# Also has a surved boundary which should be represented exactly as well 
# Except near the curved surface, the error should be on the order 
# of the number of significant digits in the .dat file (8)
# Use datatank file to examine the error.
# Not sure if curved surface should be represented exactly.  I guess not?


cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/*/*/*/*}/bin
export PATH=${PATH}:${BINDIR}

HP="tri_hp"


if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm -rf *

cp ../Inputs/* .

tri_mesh generate.inpt
mv data2_b0.grd parabola.grd
rm rstrt*
rm data*

${HP} ./run

cd ..
opendiff Results Baseline 

