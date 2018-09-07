#!/bin/bash

# This test refinement and coarsening
# with a non-uniform resolution
# Uncomment the DEBUG_ADAPT line
# in rebay.pdf to see the process for refinement
# Uncomment the DEBUG_ADAPT line
# in yaber.pdf to see the process for coarsening

cd "$(dirname "$0")"

if [ -e Results_Refine ]; then
	cd Results_Refine
else
	mkdir Results_Refine
	cd Results_Refine
fi
rm *

cp ../Inputs/* .
echo "Beginning Generate"
tri_mesh generate
echo "Ending Generate"
cd ..

opendiff Results_Refine/ Results_Refine_Baseline/
