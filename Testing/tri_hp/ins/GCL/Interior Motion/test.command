#!/bin/bash
# Tests to make sure the GCL is satisfied
# Moves mesh in x-direction and free-stream
# should stay exactly the solutin

# cd to the directory where the script resides.
# Necessary for platforms where script can be 
# double clicked from gui
cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/*/*/*/*}/bin
export PATH=${PATH}:${BINDIR}

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

# Make meshes
tri_mesh -m 'x0+0.01*x1,x1' square1.grd distorted1.grd

let ngrid=1
while [ $ngrid -lt 5 ]; do
	let ngp=${ngrid}+1
	tri_mesh -r distorted${ngrid}.grd distorted${ngp}.grd
	let ngrid=${ngrid}+1
done

let ngrid=1
while [ $ngrid -lt 6 ]; do
	tri_mesh -m 'x0-0.01*x1,x1' distorted${ngrid}.grd square${ngrid}.grd
	let ngrid=${ngrid}+1
done
rm distorted*
rm square1.grd square2.grd square3.grd square4.grd

tri_hp run.inpt

cd ..

opendiff Baseline/ Results/
