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

PARTSCRIPT="${HOME}/Codes/tet_mesh/partition.bash"

${PARTSCRIPT} run.inpt

cd ../
opendiff Baseline/ Results/
