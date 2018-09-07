#!/bin/bash
# Runs a case with a solid-fluid interface
# No flow just a change in thermal conductivity
# slope in bottom part (fluid) should be 1/3
# slope in top part.
# Mostly just a check to make sure b.c. works

# Choose this for petsc
HP="mpiexec -np 2 tri_hp_petsc run.inpt"

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
tri_mesh generate
rm data*.grd

${HP} ${PETSC}

cd ..

opendiff Baseline/ Results/
