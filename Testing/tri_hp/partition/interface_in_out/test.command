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
mod_map generate.inpt b0_mesh squarebot.d
mod_map generate.inpt b1_mesh squaretop.d
mod_map generate.inpt "growth factor" 4000
mod_map generate.inpt b0_s1_type symbolic_comm
mod_map generate.inpt b1_s1_type symbolic_comm
mod_map generate.inpt restart_interval 1
mod_map generate.inpt logfile generate
mpiexec -np 2 tri_mesh_mpi generate.inpt
rm data*.grd

#PETSC="-stop_for_debugger"
mpiexec -np 2 tri_hp_petsc run.inpt ${PETSC}

let NPART=4
mod_map run.inpt partition ${NPART}
let NTSTEP=$(mod_map -e run.inpt ntstep)
mod_map run.inpt restart $NTSTEP
tri_subpartition run.inpt ${NPART}
let TOTAL=$(mod_map -e partition.inpt nblock | wc -w | tr -d ' ')

#mpiexec -np ${TOTAL} tri_hp_petsc partition.inpt -stop_for_debugger
mpiexec -np ${TOTAL} tri_hp_petsc partition.inpt

cd ..

opendiff Baseline/ Results/
