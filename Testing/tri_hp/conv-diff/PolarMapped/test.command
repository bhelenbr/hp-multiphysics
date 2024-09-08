#!/bin/bash

# Testing accuracy for a case with a singular point
cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/*/*/*/*}/bin
export PATH=${PATH}:${BINDIR}

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm -rf *

cp ../Inputs/* .

tri_mesh generate.inpt

cp generate.inpt run.inpt
mod_map run.inpt b0_mesh rstrt1_b0.grd
mod_map run.inpt logfile run
mod_map run.inpt ncycle 20
mod_map run.inpt ntstep 1
mod_map run.inpt adapt 0


let log2p=0
while [ $log2p -lt 3 ]; do
	mkdir log2p${log2p}
	cp run.inpt log2p${log2p}
	cp generate.inpt log2p${log2p}
	cp rstrt1_b0.grd log2p${log2p}
	cd log2p${log2p}
	mod_map run.inpt log2p ${log2p}

	mpiexec -np 1 tri_hp_petsc run.inpt
	tail -2 run_b0.log | head -1 | cut -d\  -f2,4 >> cnvg.dat
		
	mod_map generate.inpt refineby2 1
	mod_map generate.inpt b0_mesh rstrt1_b0.grd
	
	let nsteps=6
	let ngrid=1
	while [ $ngrid -le $nsteps ]; do
		mod_map generate.inpt restart ${ngrid}
		mod_map generate.inpt
		mpiexec -np 1 tri_hp_petsc generate.inpt
		let ngrid=${ngrid}+1
		mod_map run.inpt restart ${ngrid}
		mpiexec -np 1 tri_hp_petsc run.inpt
		tail -2 run_b0.log | head -1 | cut -d\  -f2,4 >> cnvg.dat
		
		let ngp=${ngrid}+1
		cp rstrt${ngrid}_b0.nc rstrt${ngp}_b0.nc
	
		let ngrid=${ngrid}+1
	done
	cd ..
	let log2p=${log2p}+1
done
cd ..
./make_plot.command > Results/rates.dat
opendiff Results/ Baseline/

