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

mkdir coarsen
cd coarsen
cp ../parabola.grd .
cp ../run.inpt .
mod_map run.inpt maximum_length 0.8
mod_map run.inpt minimum_length 0.4
mod_map run.inpt logfile coarsen
${HP} ./run
cp rstrt1_d0_b0.dat ../coarsen0_b0.dat
cd ..

mkdir refine
cd refine
cp ../parabola.grd .
cp ../run.inpt .
mod_map run.inpt minimum_length 0.025
mod_map run.inpt maximum_length 0.05
mod_map run.inpt logfile refine
${HP} ./run
cp rstrt1_d0_b0.dat ../refine0_b0.dat
cd ..



cd ..
opendiff Results Baseline 

