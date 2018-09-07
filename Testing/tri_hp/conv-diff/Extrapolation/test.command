#!/bin/bash
# Test the solution extrapolation
# First two runs compare non-uniform case w/ & w/o extrapolation
# Last case should give exact answer by extrapolating

cd "$(dirname "$0")"

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm *

cp ../Inputs/* .

HP="tri_hp"

% Tests real cas
${HP} run
mod_map run.inpt extrapolate 1
${HP} run

mod_map run.inpt b0_s3_hp_type adiabatic
mod_map run.inpt b0_s4_hp_type adiabatic
mod_map run.inpt b0_s5_hp_type adiabatic
mod_map run.inpt src symbolic
mod_map run.inpt src0 1.0
${HP} run

cd ..

opendiff Baseline/ Results/
