#!/bin/bash
# Runs a case with a fixed point at the inflow
# and outflow point at the outflow

cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/*/*/*/*/*}/bin
export PATH=${PATH}:${BINDIR}

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
