#!/bin/bash
# Tests the contact angle implementation
# Solution should converge to a steady-state configuration
# with 80 degree contact angle at edges of domain

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

${HP} run

cd ..

opendiff Baseline/ Results/
