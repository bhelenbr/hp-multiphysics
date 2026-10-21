#!/bin/bash
# Plane strain: squish a unit square box (squish.inpt) and verify the
# ps module against an exact homogeneous compression solution (exact.inpt)

cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/Testing/*}/bin
export PATH=${BINDIR}:${PATH}

HP="tri_hp"

if [ -e Results ]; then
	cd Results
	rm -rf *
else
	mkdir Results
	cd Results
fi

# copy input files into results directory
cp ../Inputs/* .

# 1x1 square -> 8x8 (square4.grd)
let ngrid=1
while [ $ngrid -lt 4 ]; do
	let ngp=${ngrid}+1
	tri_mesh -r square${ngrid}.grd square${ngp}.grd
	let ngrid=${ngrid}+1
done

${HP} exact.inpt
grep L_2 exact_b0.log | tail -1

${HP} squish.inpt

cd ..
