#!/bin/bash
# Free stream flow with a sinusoidal free surface
# Runs a case with a fixed point at the inflow
# and outflow point at the outflow

cd "$(dirname "$0")"

if [ -e Results_petsc ]; then
	cd Results_petsc
else
	mkdir Results_petsc
	cd Results_petsc
fi
rm *

cp ../Inputs/* .

HP="mpiexec -np 1 tri_hp_petsc"

#PETSC="-stop_for_debugger"

${HP} run_petsc.inpt ${PETSC}

cd ..

opendiff Baseline_petsc/ Results_petsc/
