#!/bin/bash

# Generates a mesh for a curved bow shock
# where there are two domains before and after shock
#
# This actually just calculates a M= 0.3 flow
# with a communication boundary for now
# Uses cns which is based on primitive variables
# p, u, v, RT
# variables are non dimensional
# so p = 1/gamma
# u, v, = u/c
# RT = 1/gamma
# This makes rho = 1
# and c = 1

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

cd ..

#opendiff Results Baseline
