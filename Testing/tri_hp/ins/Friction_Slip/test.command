#!/bin/bash
# Runs a case with a friction slip boundary
# Domain is horizontal with symmetry on bottom
# and friction slip on top
# Like a channel flow except with slip B.C.

cd "$(dirname "$0")"

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm *

cp ../Inputs/* .

mpiexec -np 1 tri_hp_petsc run.inpt

cd ..

opendiff Baseline/ Results/
