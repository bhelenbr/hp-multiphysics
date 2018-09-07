#!/bin/bash
# Free stream flow with a sinusoidal free surface
# Runs a case with a fixed point at the inflow
# and outflow point at the outflow

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

${HP} run_mgrid

cd ..

opendiff Baseline_mgrid/ Results_mgrid/
