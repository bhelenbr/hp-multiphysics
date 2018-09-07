#!/bin/bash

# Description
# This test case matches boundaries using a single pass with edges
# Then a pass with vertices

cd "$(dirname "$0")"

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm *

cp ../Inputs/* .

tri_mesh run

cd ..
opendiff BASELINE/ Results/
