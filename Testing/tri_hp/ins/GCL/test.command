#!/bin/bash
# Tests to make sure the GCL is satisfied
# Moves mesh in x-direction and free-stream
# should stay exactly the solutin

# cd to the directory where the script resides.
# Necessary for platforms where script can be 
# double clicked from gui
cd "$(dirname "$0")"

# Make Results directory
if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm *

# copy input files into results directory
cp ../Inputs/* .

tri_hp run.inpt

cd ..

opendiff Baseline/ Results/
