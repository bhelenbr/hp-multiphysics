#!/bin/bash
# Translation test case where solution translates
# and rigid boundary translates as well.(translating_surface b.c.) 
# plot results to see that everything is in sync.

cd "$(dirname "$0")"

if [ -e Results ]; then
	rm -rf Results
fi
mkdir Results
cd Results

cp ../Inputs/* .

mpiexec -np 1 tri_hp_petsc run.inpt

cd ..
opendiff Results/ Baseline/
