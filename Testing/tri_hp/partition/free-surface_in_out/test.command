#!/bin/bash
# Runs a case with a fixed point at the inflow
# and outflow point at the outflow
# I don't think jacobian at outflow point is correct

cd "$(dirname "$0")"

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm *

cp ../Inputs/* .

cp run.inpt generate.inpt
mod_map generate.inpt ntstep 1
mod_map generate.inpt adapt 1
mod_map generate.inpt b0_mesh square.d
mod_map generate.inpt "growth factor" 4000
mod_map generate.inpt b0_s3_type symbolic
mod_map generate.inpt restart_interval 1
mod_map generate.inpt logfile generate

tri_mesh generate.inpt
rm data*.grd

#PETSC="-stop_for_debugger"
mpiexec -np 1 tri_hp_petsc run.inpt ${PETSC}

let NPART=4
mod_map run.inpt partition ${NPART}
let NTSTEP=$(mod_map -e run.inpt ntstep)
mod_map run.inpt restart $NTSTEP
~/Codes/tri_hp/tri_subpartition run.inpt ${NPART}

mpiexec -np ${NPART} tri_hp_petsc partition.inpt

cd ..

opendiff Baseline_petsc/ Results/
