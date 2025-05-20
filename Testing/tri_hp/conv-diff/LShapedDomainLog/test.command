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

# set -e  
tri_mesh generate.inpt


# cp generate.inpt run.inpt
# mod_map run.inpt b0_mesh rstrt1_b0.grd
# mod_map run.inpt b1_mesh rstrt1_b1.grd
# # mod_map run.inpt "growth factor" 10
# mod_map run.inpt b0_type cd
# mod_map run.inpt b1_type cd
# mod_map run.inpt nblock "1 1"
# mod_map run.inpt logfile run
# mod_map run.inpt adapt 0
# mod_map run.inpt ntstep 1
# mod_map run.inpt ncycle 20




# let log2p=0

# mkdir log2p${log2p}
# cp run.inpt log2p${log2p}
# cp generate.inpt log2p${log2p}
# cp rstrt1_b0.grd log2p${log2p}
# cp rstrt1_b1.grd log2p${log2p}
# cd log2p${log2p}
# mod_map run.inpt log2p ${log2p}

# mpiexec -np 2 tri_hp_petsc run.inpt
# tail -2 run_b0.log | head -1 | cut -d\  -f2,4 | tr -d '\n' >> cnvg.dat
# echo -n ' ' >> cnvg.dat
# grep DOF run_b0.log | cut -d\  -f6 >> cnvg.dat
# mod_map generate.inpt refineby2 1
# mod_map generate.inpt b0_mesh rstrt1_b0.grd
# mod_map generate.inpt b1_mesh rstrt1_b1.grd

# let nsteps=6
# let ngrid=1

# mod_map generate.inpt restart ${ngrid}
# mod_map generate.inpt
# mpiexec -np 2 tri_hp_petsc generate.inpt




# let ngrid=${ngrid}+1
# mod_map run.inpt restart ${ngrid}
# mpiexec -np 2 tri_hp_petsc run.inpt
# tail -2 run_b0.log | head -1 | cut -d\  -f2,4 | tr -d '\n' >> cnvg.dat
# echo -n ' ' >> cnvg.dat
# grep DOF run_b0.log | cut -d\  -f6 >> cnvg.dat
# let ngp=${ngrid}+1
# cp rstrt${ngrid}_b0.nc rstrt${ngp}_b0.nc

# let ngrid=${ngrid}+1

# while [ $log2p -lt 3 ]; do
# 	mkdir log2p${log2p}
# 	cp run.inpt log2p${log2p}
# 	cp generate.inpt log2p${log2p}
# 	cp rstrt1_b0.grd log2p${log2p}
# 	cp rstrt1_b1.grd log2p${log2p}
# 	cd log2p${log2p}
# 	mod_map run.inpt log2p ${log2p}

# 	mpiexec -np 2 tri_hp_petsc run.inpt
# 	tail -2 run_b0.log | head -1 | cut -d\  -f2,4 | tr -d '\n' >> cnvg.dat
# 	echo -n ' ' >> cnvg.dat
# 	grep DOF run_b0.log | cut -d\  -f6 >> cnvg.dat
# 	mod_map generate.inpt refineby2 1
# 	mod_map generate.inpt b0_mesh rstrt1_b0.grd
# 	mod_map generate.inpt b1_mesh rstrt1_b1.grd
	
# 	let nsteps=6
# 	let ngrid=1
# 	while [ $ngrid -le $nsteps ]; do
# 		mod_map generate.inpt restart ${ngrid}
# 		mod_map generate.inpt
# 		mpiexec -np 2 tri_hp_petsc generate.inpt
# 		let ngrid=${ngrid}+1
# 		mod_map run.inpt restart ${ngrid}
# 		mpiexec -np 2 tri_hp_petsc run.inpt
# 		tail -2 run_b0.log | head -1 | cut -d\  -f2,4 | tr -d '\n' >> cnvg.dat
# 		echo -n ' ' >> cnvg.dat
# 		grep DOF run_b0.log | cut -d\  -f6 >> cnvg.dat
# 		let ngp=${ngrid}+1
# 		cp rstrt${ngrid}_b0.nc rstrt${ngp}_b0.nc
	
# 		let ngrid=${ngrid}+1
# 	done
# 	cd ..
# 	let log2p=${log2p}+1
# done
# cd ..
# ./make_plot.command > Results/rates.dat
# opendiff Results/ Baseline/



# let log2p=0
# while [ $log2p -lt 3 ]; do
# 	mkdir log2p${log2p}
# 	cp run.inpt log2p${log2p}
# 	cp offset.inpt log2p${log2p}
# 	cp rstrt3_b0.grd log2p${log2p}
# 	cp rstrt3_b1.grd log2p${log2p}
# 	cd log2p${log2p}
# 	mod_map run.inpt log2p ${log2p}

# 	mpiexec -np 2 tri_hp_petsc run.inpt
# 	tail -2 output_b0.log | head -1 | cut -d\  -f2,4 | tr -d '\n' >> cnvg.dat
# 	echo -n ' ' >> cnvg.dat
# 	grep DOF output_b0.log | cut -d\  -f6 >> cnvg.dat
# 	mod_map offset.inpt refineby2 1
# 	mod_map offset.inpt b0_mesh rstrt3_b0.grd
# 	mod_map offset.inpt b1_mesh rstrt3_b1.grd
	
# 	let nsteps=6
# 	let ngrid=1
# 	while [ $ngrid -le $nsteps ]; do
# 		mod_map run.inpt restart ${ngrid}
# 		mpiexec -np 2 tri_hp_petsc offset.inpt
# 		let ngrid=${ngrid}+1
# 		mod_map run.inpt restart ${ngrid}
# 		mpiexec -np 2 tri_hp_petsc run.inpt
# 		tail -2 output_b0.log | head -1 | cut -d\  -f2,4 | tr -d '\n' >> cnvg.dat
# 		echo -n ' ' >> cnvg.dat
# 		grep DOF output_b0.log | cut -d\  -f6 >> cnvg.dat
# 		let ngp=${ngrid}+1
# 		cp rstrt${ngrid}_b0.nc rstrt${ngp}_b0.nc
	
# 		let ngrid=${ngrid}+1
# 	done
# 	cd ..
# 	let log2p=${log2p}+1
# done
# cd ..
# ./make_plot.command > Results/rates.dat
# opendiff Results/ Baseline/

