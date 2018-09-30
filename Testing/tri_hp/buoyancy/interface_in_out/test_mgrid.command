#!/bin/bash
# Tests the same problem as ins/interface_in_out
# except makes sure that this works when we add the
# temperature as a continuous unknown

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

tri_hp run_mgrid.inpt

cd ..

opendiff Baseline_mgrid/ Results_mgrid/
