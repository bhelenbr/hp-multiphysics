#!/bin/bash
# Tests the same problem as ins/interface_in_out
# except makes sure that this works when we add the
# temperature as a continuous unknown

cd "$(dirname "$0")"

if [ -e Results_mgrid ]; then
	cd Results_mgrid
else
	mkdir Results_mgrid
	cd Results_mgrid
fi
rm *

cp ../Inputs/* .

HP="$HOME/Codes/tri_hp/build/Release/tri_hp"
HP="$HOME/bin/tri_hp"

${HP} run_mgrid.inpt

cd ..

opendiff Baseline_mgrid/ Results_mgrid/
