#!/bin/bash
# Tests the same problem as ins/interface_in_out
# except makes sure that this works when we add the
# temperature as a continuous unknown

cd "$(dirname "$0")"

if [ -e Results_petsc ]; then
	cd Results_petsc
else
	mkdir Results_petsc
	cd Results_petsc
fi
rm *

cp ../Inputs/* .

HP="mpiexec -np 2 tri_hp_petsc"

#PETSC="-stop_for_debugger"

${HP} run_petsc.inpt ${PETSC}

cd ..

opendiff Baseline_petsc/ Results_petsc/
