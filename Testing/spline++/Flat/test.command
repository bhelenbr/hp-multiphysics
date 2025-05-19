#!/bin/bash
# This is a test of spline interpolation for a straight line

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

# generate mesh and remove unnecessary data files
spline -i pts.dat flat.spl > results.dat
spline -o 0.1 -i pts.dat flat.spl >> results.dat

cd ..

# use opendiff (on OS X) to compare to a Baseline run
# can change to diff to do this on linux
opendiff Baseline/ Results/
