#!/bin/bash
# Runs a free-surface case with periodic boundary conditions
# Check convergence and check data files to make sure 
# periodic points have exactly the same values

cd "$(dirname "$0")"

if [ -e Results_mgrid ]; then
	cd Results_mgrid
else
	mkdir Results_mgrid
	cd Results_mgrid
fi
rm *

cp ../Inputs/* .

HP="tri_hp"

${HP} run_mgrid.inpt

cd ..

opendiff Baseline_mgrid/ Results_mgrid/
