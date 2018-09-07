#!/bin/bash

# This tests the phased communication 
# Two sweeps 1 horizontal then 1 vertical to get vertex right

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
