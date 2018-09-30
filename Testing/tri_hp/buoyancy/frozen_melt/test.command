#!/bin/bash
# Runs a case with a solid-fluid interface where the 
# temperature is fixed on the interface.  This is to
# set-up initial conditions for a melting problem

# cd to the directory where the script resides.
# Necessary for platforms where script can be 
# double clicked from gui
cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/*/*/*/*}/bin
export PATH=${PATH}:${BINDIR}

HP="mpiexec -np 2 tri_hp_petsc run.inpt"

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
