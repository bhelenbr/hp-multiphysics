#!/bin/bash

# Generates a mesh for a nozzle
#
# Solves steady state problem, can change back pressure in 
# cns/getnewibc.cpp to make shock move
# 
# Note that downstream M in run.inpt comes from running the nozzle
# ibc with meaningless numbers in the downstream and getting
#  the resulting M value at the shock location

# This runs a nozzle with shock at xs = 7. Uses a modified version of the nozzle shape in https://www.grc.nasa.gov/www/wind/valid/cdv/cdv.html

# Nozzle_fun.m is used to determine the nozzle shape parameters.
# Nozzle_shape.m plots nozzle shape and calculates nozzle height at an x location.
# Nozzle_yvel.m outlines how the y-velocity is initialized in the nozzle.

# To run with a different configuration (different xs, etc) you may need to change M in run.inpt #IC Downstream.
# This is the Mach number on the upstream side of the shock. Get this by running the test case with any number
# greater than 1 for M, and output the correct M from Compressible/getnewibc::Nozzle Newton iteration loop.
# Then re-run with the correct value of M on the upstream side of the shock.

# This only seems to work with stabilization off

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

opendiff Results Baseline
