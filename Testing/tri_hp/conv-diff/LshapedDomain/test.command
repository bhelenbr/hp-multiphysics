#!/bin/bash

# Testing accuracy for a case with a singular point
cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/Testing/*}/bin
echo ${BINDIR}
export PATH=${BINDIR}:${PATH}

set -e

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
mod_map run.inpt b0_mesh rstrt0_b0.grd
mod_map run.inpt b1_mesh rstrt0_b1.grd
mod_map run.inpt "growth factor" 10
mod_map run.inpt b0_type cd
mod_map run.inpt b1_type cd
mod_map run.inpt nblock "1 1"
mod_map run.inpt logfile output
mod_map run.inpt adapt 0
mod_map run.inpt ncycle 10
mod_map run.inpt ntstep 1

log2p=0

while [ $log2p -lt 3 ]; do
	mkdir log2p${log2p}
	cp run.inpt log2p${log2p}
	cp rstrt3_b0.grd log2p${log2p}/rstrt0_b0.grd
	cp rstrt3_b1.grd log2p${log2p}/rstrt0_b1.grd
	cd log2p${log2p}
	mod_map run.inpt log2p ${log2p}

	mpiexec -np 2 tri_hp_petsc run.inpt
	
	dof0=$(grep 'DOF:' output_b0.log | cut -d' ' -f6)
	dof1=$(grep 'DOF:' output_b1.log | cut -d' ' -f6)
	((DOF = dof0 + dof1))
	echo ${DOF} | tr -d '\n' >> cnvg.dat
	echo -n ' ' >> cnvg.dat
	tail -1 output_b0.log | cut -d\  -f3 | tr -d '\n' >> cnvg.dat
	echo -n ' ' >> cnvg.dat
	tail -2 output_b0.log | head -1 | cut -d\  -f2,4 >> cnvg.dat

	grep '#L_2' output_b0.log | head -1 | cut -d\  -f2,4 >> ic.dat 

	ngrids=5
	ngrid=1
	restart=1
	while [ $ngrid -lt $ngrids ]; do
		# Refine solution
		mod_map run.inpt refineby2 1
		mod_map run.inpt adapt 1
		mod_map run.inpt restart ${restart}
		mpiexec -np 2 tri_hp_petsc run.inpt
		((restart++))
		
		# Run case using restart file
		mod_map run.inpt restart ${restart}
		mod_map run.inpt adapt 0
		mod_map run.inpt refineby2 0
		mod_map run.inpt b0_mesh rstrt${restart}_b0.nc
		mod_map run.inpt b1_mesh rstrt${restart}_b1.nc
		mpiexec -np 2 tri_hp_petsc run.inpt
		((restart++))
		
		let DOF=$(grep 'DOF:' output_b0.log | cut -d\  -f6)+$(grep 'DOF:' output_b1.log | cut -d\  -f6)
		dof0=$(grep 'DOF:' output_b0.log | cut -d' ' -f6)
		dof1=$(grep 'DOF:' output_b1.log | cut -d' ' -f6)
		((DOF = dof0 + dof1))
		echo ${DOF} | tr -d '\n' >> cnvg.dat
		echo -n ' ' >> cnvg.dat
		tail -1 output_b0.log | cut -d\  -f3 | tr -d '\n' >> cnvg.dat
		echo -n ' ' >> cnvg.dat
		tail -2 output_b0.log | head -1 | cut -d\  -f2,4 >> cnvg.dat
		
		((ngrid++))
	done
	
	mkdir ICtest
	cd ICtest
	cp ../run.inpt .
	cp ../rstrt*_b?.nc .
	# Test initial conditions
	ngrid=1
	restart=0
	mod_map -d run.inpt restart
	mod_map run.inpt ncycle 0
	while [ $ngrid -lt $ngrids ]; do
		((restart += 2))
		
		# Run case using restart file
		mod_map run.inpt b0_mesh rstrt${restart}_b0.nc
		mod_map run.inpt b1_mesh rstrt${restart}_b1.nc
		mpiexec -np 2 tri_hp_petsc run.inpt
		grep '#L_2' output_b0.log | head -1 | cut -d\  -f2,4 >> ../ic.dat 
		mv data0_b0.dat data${restart}_b0.dat
		mv data0_b1.dat data${restart}_b1.dat
		rm data1_b0.dat
		rm data1_b1.dat
		((ngrid++))
	done
	cd ../..
	((log2p++))
done
cd ..
./make_plot.command > Results/rates.dat
opendiff Results/ Baseline/

