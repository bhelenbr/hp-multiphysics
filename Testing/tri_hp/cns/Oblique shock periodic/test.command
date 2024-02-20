#!/bin/bash

# Runs a quasi-1D shock problem
# Shock is a vertical line separating upstream (left)
# and downstream (right) blocks of the domain

cd "$(dirname "$0")"
BINDIR=${PWD%/*/*/*/*}/bin
export PATH=${PATH}:${BINDIR}

HP="tri_hp_petsc"
#PETSC="-stop_for_debugger"

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm *

cp ../Inputs/* .

# Create structured grid file


mpiexec -np 2 ${HP} run.inpt ${PETSC}

cd ..

opendiff Results Baseline
