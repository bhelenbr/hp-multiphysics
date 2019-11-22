#!/bin/bash

# Generates a mesh for a nozzle
#
# Solves steady state problem, can change back pressure in 
# cns/getnewibc.cpp to make shock move
# 
# Note that downstream M in run.inpt comes from running the nozzle
# ibc with meaningless numbers in the downstream and getting
#  the resulting M value at the shock location

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

tri_mesh generate

mpiexec -np 2 ${HP} run.inpt ${PETSC}

# mod_map run.inpt adapt 1
# mod_map run.inpt error_estimator energy_norm
# mod_map run.inpt error_target 1e-4
# mod_map run.inpt adapt_output 1
# 
# mpiexec -np 2 ${HP} run.inpt ${PETSC}

#cd ..

#opendiff Results Baseline
