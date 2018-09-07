#!/bin/bash

# This tests refinement and coarsening
# with a non-uniform resolution
# Uncomment the DEBUG_ADAPT line
# in rebay.pdf to see the process for refinement
# Uncomment the DEBUG_ADAPT line
# in yaber.pdf to see the process for coarsening

cd "$(dirname "$0")"

if [ -e Results_Coarsen ]; then
	cd Results_Coarsen
else
	mkdir Results_Coarsen
	cd Results_Coarsen
fi
rm *

cp ../Inputs/* .

# MAKE A NEW LENGTH INPUT FILE SO THAT tri_mesh GETS COARSENED
NVRTX0=$(head -1 generate_b0.grd | cut -d" " -f 2)
COUNTER=0
while [  $COUNTER -lt $NVRTX0 ]; do
 echo "0.25" >> generate_b0.lngth
 let COUNTER=COUNTER+1 
done

echo "Beginning Coarsen"
tri_mesh coarsen
echo "Ending Coarsen"

cd ..

opendiff Results_Coarsen/ Results_Coarsen_Baseline/
