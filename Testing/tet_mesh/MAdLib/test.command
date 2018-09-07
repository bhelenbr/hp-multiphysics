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

MSH="$HOME/bin/tet_mesh"

${MSH} run

cd ../
opendiff Baseline/ Results/
