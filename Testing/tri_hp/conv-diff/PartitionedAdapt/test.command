#!/bin/bash
# Runs an adapting case in parallel with known exact solution & geometry
# No iteration neessay
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

tri_mesh generate.inpt
mv rstrt1_b0.grd parabola.grd
rm data*.grd

mkdir coarsen
cd coarsen
cp ../run.inpt .
cp ../parabola.grd .
mod_map run.inpt maximum_length 0.8
mod_map run.inpt minimum_length 0.4
mpiexec -np 1 tri_hp_petsc run.inpt

let NPART=4
mod_map run.inpt partition ${NPART}
mod_map run.inpt restart 2
mod_map run.inpt ntstep 1
mod_map run.inpt logfile partition_log
tri_subpartition run.inpt ${NPART}
if [ "$?" -ne "0" ]; then
        echo "partitioning failed"
        exit 1
fi
mod_map partition.inpt ntstep 2
mpiexec -np ${NPART} tri_hp_petsc partition.inpt ${PETSC}
cd ..

exit 1

mkdir refine
cd refine
cp ../run.inpt .
cp ../parabola.grd .
mod_map run.inpt minimum_length 0.025
mod_map run.inpt maximum_length 0.05
mpiexec -np 1 tri_hp_petsc run.inpt

let NPART=4
mod_map run.inpt partition ${NPART}
mod_map run.inpt restart 2
mod_map run.inpt ntstep 1
mod_map run.inpt logfile partition_log
tri_subpartition run.inpt ${NPART}
if [ "$?" -ne "0" ]; then
        echo "partitioning failed"
        exit 1
fi
mod_map partition.inpt ntstep 2
mpiexec -np ${NPART} tri_hp_petsc partition.inpt ${PETSC}
cd ..

cd ..

opendiff Baseline/ Results/
