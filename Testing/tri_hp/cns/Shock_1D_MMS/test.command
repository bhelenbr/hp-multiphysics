#!/bin/bash

# This runs a vertical fully developed channel (periodic flow driven by gravity)

cd "$(dirname "$0")"
BINDIR=${PWD%/*/*/*/*}/bin
export PATH=${PATH}:${BINDIR}
HP="tri_hp_petsc_MMS"

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm *

#PETSC="-stop_for_debugger"

cp ../Inputs/* .

tri_mesh generate
mv rstrt1_b0.grd square1_b0.grd
mv rstrt1_b1.grd square1_b1.grd

let nrefinements=3
let log2p=2
let tstep=200
#while [ $log2p -lt 3 ]; do
	let ngrid=1
	while [ $ngrid -le $nrefinements ]; do
		mod_map run.inpt b0_mesh square${ngrid}_b0.grd
		mod_map run.inpt b1_mesh square${ngrid}_b1.grd
		mod_map run.inpt log2p ${log2p}
		mod_map run.inpt logfile output${ngrid}_p${log2p}
		mod_map run.inpt dtinv cu_avg*10*${tstep}/_pi
		mod_map run.inpt ntstep 0.001*dtinv*2*_pi/omegas
		mpiexec -np 2 ${HP} run.inpt ${PETSC}
		tail -2 output${ngrid}_p${log2p}_b0.log | head -1 | cut -d\  -f2,4,6,8,10,12,14,16 >> cnvg${log2p}_b0.dat
		tail -2 output${ngrid}_p${log2p}_b1.log | head -1 | cut -d\  -f2,4,6,8,10,12,14,16 >> cnvg${log2p}_b1.dat
		let ngp=${ngrid}+1
		tri_mesh -r square${ngrid}_b0.grd square${ngp}_b0.grd
		tri_mesh -r square${ngrid}_b1.grd square${ngp}_b1.grd
		let ngrid=${ngp}
		let tstep=tstep*8
	done
	let log2p=${log2p}+1
	let tstep=200
#done

cd ..

#opendiff Results Baseline
