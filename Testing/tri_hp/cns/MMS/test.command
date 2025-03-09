#!/bin/bash

# This runs The method of manufactured solutions for compressible Navier Stokes with 
# properties sinusoidally oscillating in space

# To run uncomment Define MMS in tri_hp_cns.h


cd "$(dirname "$0")"
BINDIR=${PWD%/*/*/*/*}/bin
export PATH=${PATH}:${BINDIR}
HP="tri_hp_petsc"
#PETSC="-stop_for_debugger"

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm *

cp ../Inputs/* .

tri_mesh generate
tri_mesh -r rstrt1_b0.grd square1_b0.grd
mv rstrt1_b0.grd square1_b0.grd

let nrefinements=3
let log2p=0
while [ $log2p -lt 3 ]; do
	let ngrid=1
	while [ $ngrid -le $nrefinements ]; do
		mod_map run.inpt b0_mesh square${ngrid}_b0.grd
		mod_map run.inpt log2p ${log2p}
		mod_map run.inpt logfile output${ngrid}_p${log2p}
		mpiexec -np 1 ${HP} run.inpt ${PETSC}
		tail -2 output${ngrid}_p${log2p}_b0.log | head -1 | cut -d\  -f2,4,6,8,10,12,14,16 >> cnvg${log2p}.dat
		let ngp=${ngrid}+1
		tri_mesh -r square${ngrid}_b0.grd square${ngp}_b0.grd
		let ngrid=${ngp}
	done
	let log2p=${log2p}+1
done

cd ..

./make_plot.command > Results/rates.dat

opendiff Results Baseline
