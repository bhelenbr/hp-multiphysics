#!/bin/bash

## Test of partitioning script

# cd to the directory where the script resides.
# Necessary for platforms where script can be 
# double clicked from gui
cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/*/*/*}/bin
echo ${BINDIR}
export PATH=${PATH}:${BINDIR}

NPROC=(2 4 8)

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm -rf *

cp ../Inputs/* .
tri_mesh -A rstrt1_b0.grd rstrt1_b1.grd merge.grd

let npc=0
while [ $npc -lt ${#NPROC[@]} ]; do
	let np=${NPROC[npc]}
	mkdir partition${np}
	cd partition${np}
	cp ../* .
	tri_mesh -p merge.grd ${NPROC[$npc]}
	cd ..
	let npc=$npc+1
done

cd ..

opendiff Baseline/ Results/
