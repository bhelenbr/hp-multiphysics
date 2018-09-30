#!/bin/bash
# Runs a free-surface case with periodic boundary conditions
# Check convergence and check data files to make sure 
# periodic points have exactly the same values
# This should work whether or not the vertex points
# are defined as periodic

cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/*/*/*/*/*}/bin
export PATH=${PATH}:${BINDIR}

if [ -e Results_petsc ]; then
	cd Results_petsc
else
	mkdir Results_petsc
	cd Results_petsc
fi
rm *

cp ../Inputs/* .

#HP="$HOME/Desktop/tri_hp1/build/Release/tri_hp"

HP="mpiexec -np 1 tri_hp_petsc"

${HP} run_petsc.inpt

cd ..

opendiff Baseline_petsc/ Results_petsc/
