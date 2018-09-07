#!/bin/bash
# Coarsens a mesh
# For some reason this case doesn't work MAdLib crashes
# At one point it worked if you output the mesh and re-input the mesh
# lines must be uncommented in 
# MAdLibInterface::coarsenMesh
# to make it work (still doesn't)

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

tet_mesh run

cd ../
opendiff Baseline/ Results/
