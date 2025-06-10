#!/bin/bash

# Testing consistency for a spline mapping
# Should give 0 for L2 error

cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/Testing/*}/bin
export PATH=${BINDIR}:${PATH}

set -e

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm -rf *

cp ../Inputs/* .

tri_mesh generate.inpt


cp generate.inpt run.inpt
mod_map run.inpt b0_mesh rstrt2_b0.grd
mod_map run.inpt logfile run
mod_map run.inpt ncycle 20
mod_map run.inpt ntstep 1
mod_map run.inpt adapt 0

mpiexec -np 1 tri_hp_petsc run.inpt

cd ..

opendiff Results/ Baseline/