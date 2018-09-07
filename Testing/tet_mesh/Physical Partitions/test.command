#!/bin/bash

cd "$(dirname "$0")"

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm *

# Solve for the modes
cp ../Inputs/* .

tet_mesh -gp sharedface.msh 2 > junk.dat

cd ../
opendiff Baseline/ Results/
