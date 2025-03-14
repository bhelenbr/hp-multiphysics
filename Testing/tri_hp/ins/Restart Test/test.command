#!/bin/bash

# This tests the restart capabilities to make sure
# restarting doesn't cause changes in the solution
# binio stopped working and I don't know why
# I think this is why I switched to netcdf
# check diff files - should be empty for success

cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/*/*/*/*}/bin
export PATH=${PATH}:${BINDIR}

HP="tri_hp"
export HP

P="1 2 4"
#P="1"
MGFLAG=0
#TYPES=(1 2 7)
TYPES=(1 7)

export P
export MGFLAG

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm -rf *

cp ../Inputs/*.grd .
# Make meshes
tri_mesh -m 'x0+0.01*x1,x1' square1.grd distorted1.grd

let ngrid=1
while [ $ngrid -lt 3 ]; do
	let ngp=${ngrid}+1
	tri_mesh -r distorted${ngrid}.grd distorted${ngp}.grd
	let ngrid=${ngrid}+1
done

let ngrid=1
while [ $ngrid -lt 4 ]; do
	tri_mesh -m 'x0-0.01*x1,x1' distorted${ngrid}.grd square${ngrid}.grd
	let ngrid=${ngrid}+1
done
rm distorted*


let n=0
while [ ${n} -lt ${#TYPES[@]} ]; do
	# Test just flow mesh restarts
	echo Flow${TYPES[n]}
	mkdir Flow${TYPES[n]}
	cd Flow${TYPES[n]}
	cp ../../Inputs/* .
	# 1,2,6 is txt,bin,nc
	mod_map run.inpt restart_type ${TYPES[n]}
	mod_map run.inpt reload_type ${TYPES[n]}
	./tests.bash
	cd ../

	# Test deforming mesh restarts
	echo Deforming${TYPES[n]}
	mkdir Deforming${TYPES[n]}
	cd Deforming${TYPES[n]}
	cp ../../Inputs/* .
	mod_map -u run.inpt mesh_movement
	mod_map run.inpt restart_type ${TYPES[n]}
	mod_map run.inpt reload_type ${TYPES[n]}
	./tests.bash
	cd ../
	
	# Test adapting mesh restarts
	echo Adapting${TYPES[n]}
	mkdir Adapting${TYPES[n]}
	cd Adapting${TYPES[n]}
	cp ../../Inputs/* .
	mod_map -u run.inpt mesh_movement
	mod_map run.inpt adapt 1 
	mod_map run.inpt restart_type ${TYPES[n]}
	mod_map run.inpt reload_type ${TYPES[n]}
	./tests.bash
	cd ../
	
	# Test curved mesh restarts
	echo Curved${TYPES[n]}
	mkdir Curved${TYPES[n]}
	cd Curved${TYPES[n]}
	cp ../../Inputs/* .
	mod_map -u run.inpt mesh_movement
	mod_map run.inpt adapt 1 
	mod_map run.inpt restart_type ${TYPES[n]}
	mod_map run.inpt reload_type ${TYPES[n]}
	mod_map -u run.inpt b0_s4_type
	mod_map -u run.inpt b0_s4_curved
	cp run.inpt square_bdry.inpt
	cp ../square2.grd square.grd
	tri_mesh -r square curved.grd
	tri_mesh -f curved curved
	mod_map run.inpt b0_mesh ../curved
	./tests.bash
	cd ../

	let n=$n+1
done