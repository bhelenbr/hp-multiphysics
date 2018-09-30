#!/bin/bash
# Runs a case with a solidifying surface
# Interface is planar, analytic solution is known
# This case is not working.  See melt_with_solid for info.
# Baseline_mgrid is very old

# To run_petsc test case comment 
NORMAL=true

cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/*/*/*/*}/bin
export PATH=${PATH}:${BINDIR}

if [ -e Results_mgrid ]; then
	cd Results_mgrid
else
	mkdir Results_mgrid
	cd Results_mgrid
fi
rm *

cp ../Inputs/* .

tri_mesh generate
rm data*.grd

HP="mpiexec -np 2 tri_hp_mpi run_mgrid.inpt ${PETSC}"

if [ -n "$NORMAL" ]; then
	if [ -e NORMAL ]; then
		cd NORMAL
	else
		mkdir NORMAL
		cd NORMAL
	fi
	rm *

	cp ../*.grd .
	cp ../run_mgrid.inpt .
	${HP}
	cd ..
fi
	
cd ..

opendiff Baseline_mgrid/ Results_mgrid/
