#!/bin/bash
# Runs an adapting case in parallel with free-surface movement
# and partitioning

cd "$(dirname "$0")"

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm *

cp ../Inputs/* .

VALGRIND_FLAGS+=" --track-origins=yes"
VALGRIND_FLAGS+=" --leak-check=full"
VALGRIND_FLAGS+=" --dsymutil=yes"

cp run.inpt generate.inpt
mod_map generate.inpt ntstep 1
mod_map generate.inpt ncycle 0
mod_map generate.inpt adapt 1
mod_map generate.inpt b0_mesh square.d
mod_map generate.inpt "growth factor" 4000
mod_map generate.inpt b0_s3_type symbolic
mod_map generate.inpt restart_interval 1
mod_map generate.inpt logfile generate

tri_mesh generate.inpt
rm data*.grd

mod_map run.inpt adapt 0
mod_map run.inpt ntstep 1
#PETSC="-stop_for_debugger"
mpiexec -np 1 tri_hp_petsc run.inpt ${PETSC}
#exit 1
#mod_map run.inpt ntstep 1
#mod_map run.inpt restart 1
#mpiexec -np 1 tri_hp_petsc run.inpt ${PETSC}
#mpiexec -np 1 valgrind ${VALGRIND_FLAGS} tri_hp_petsc run.inpt 
#exit


let NPART=4
mod_map run.inpt partition ${NPART}
let NTSTEP=$(mod_map -e run.inpt ntstep)
mod_map run.inpt restart $NTSTEP
mod_map run.inpt ntstep 1
mod_map run.inpt adapt 1
tri_subpartition run.inpt ${NPART}
if [ "$?" -ne "0" ]; then
        echo "partitioning failed"
        exit 1
fi

#PETSC="-stop_for_debugger"
mpiexec -np ${NPART} tri_hp_petsc partition.inpt ${PETSC}
#mpiexec -np ${NPART} valgrind ${VALGRIND_FLAGS} tri_hp_petsc partition.inpt 

cd ..

opendiff Baseline/ Results/
