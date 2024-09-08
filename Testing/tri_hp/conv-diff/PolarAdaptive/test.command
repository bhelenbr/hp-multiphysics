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
mod_map run.inpt ntstep 10


let log2p=1
while [ $log2p -lt 3 ]; do
	mkdir log2p${log2p}
	cp run.inpt log2p${log2p}
	cp generate.inpt log2p${log2p}
	cp rstrt1_b0.grd log2p${log2p}
	cd log2p${log2p}
	mod_map run.inpt log2p ${log2p}

	mpiexec -np 1 tri_hp_petsc run.inpt
	tail -20 run_b0.log | grep L_2 | cut -d\  -f2,4 | tr -d '\n' >> cnvg.dat
	echo -n ' ' >> cnvg.dat
	tail -20 run_b0.log | grep '# DOF:' | cut -d\  -f3,5,8,10 >> cnvg.dat
	let nsteps=5+$log2p
	let ngrid=2
	let nstart=10
	while [ $ngrid -le $nsteps ]; do
		mod_map run.inpt n ${ngrid}
		mod_map run.inpt restart ${nstart}
		mpiexec -np 1 tri_hp_petsc run.inpt
		tail -20 run_b0.log | grep L_2 | cut -d\  -f2,4 | tr -d '\n' >> cnvg.dat
		echo -n ' ' >> cnvg.dat
		tail -20 run_b0.log | grep '# DOF:' | cut -d\  -f3,5,8,10 >> cnvg.dat
		let ngrid=${ngrid}+1
		let nstart=${nstart}+10
	done
	cd ..
	let log2p=${log2p}+1
done
cd ..
./make_plot.command > Results/rates.dat
opendiff Results/ Baseline/

