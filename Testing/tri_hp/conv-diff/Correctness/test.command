#!/bin/bash

# Solves diffusion problem with quadratic elements with source term
# tests source term implementation and diffusion terms
# result should be quadratic function
# With curved off and log2p = 1 should get exact answer
# With curved on and log2p = 1 will have to iterate

cd "$(dirname "$0")"

# Make Results directory
if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm -rf *

cp ../Inputs/* .


HP="tri_hp"


${HP} run.inpt

mod_map run.inpt ax 1.0
mod_map run.inpt nu 0.0
mod_map run.inpt ibc0 "2*(x0-8)"
${HP} run.inpt

cd ..

opendiff Baseline/ Results/
