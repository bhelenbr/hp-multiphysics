#!/bin/bash
#  Calculates a flat plate boundary and tests partitioning and parallel speedup

cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/*/*/*/*}/bin
export PATH=${PATH}:${BINDIR}
HP="tri_hp_petsc"

rm -rf Results
mkdir Results
cd Results
cp ../Inputs/* .

# Make basic mesh
tri_mesh generate.inpt

#PETSC="-stop_for_debugger"

OMPIFLAGS="--bind-to socket"
#MPICHFLAGS="-bind-to rr"
#MPIFLAGS=${OMPIFLAGS}
MPIFLAGS=${MPICHFLAGS}
export OMP_NUM_THREADS=1
export OMP_DISPLAY_AFFINITY=true


let NREFINE=1
let NREFINEMAX=2
while [ ${NREFINE} -le ${NREFINEMAX} ]; do
	mkdir nrefine${NREFINE}
	cd nrefine${NREFINE}
	cp ../* .
	
	# RUN SINGLE PROCESS CASE
	CMD="mpiexec -np 1 ${MPIFLAGS} ${HP} run.inpt"
	echo ${CMD}
	eval "${CMD}"
	echo -n "1 " >> cputimes.dat
	tail -1 out_b0.log | cut -d \  -f 3 | tr -d '\n' >> cputimes.dat
	echo -n " " >> cputimes.dat
	grep 'jacobian made' out_b0.log | tail -1 | cut -d\  -f 3 | tr -d '\n' >> cputimes.dat 
	echo -n " " >> cputimes.dat
	grep 'matrix inverted' out_b0.log | tail -1 | cut -d\  -f 3 >> cputimes.dat
	
	let NPART=2
	let NPARTMAX=8
	while [ ${NPART} -le ${NPARTMAX} ]; do
		mkdir npart${NPART}
		cd npart${NPART}
		cp ../* .
		mod_map run.inpt adapt 1
		tri_subpartition run.inpt ${NPART}
		if [ "$?" -ne "0" ]; then
				echo "partitioning failed"
				exit 1
		fi
		mod_map partition.inpt adapt 0
		CMD="mpiexec -np ${NPART} ${MPIFLAGS} ${HP} partition.inpt"
		echo ${CMD}
		eval "${CMD}"
		echo -n "${NPART} "  >> ../cputimes.dat
		tail -1 out_b0.log | cut -d \  -f 3 | tr -d '\n' >> ../cputimes.dat
		echo -n " " >> ../cputimes.dat
		grep 'jacobian made' out_b0.log | tail -1 | cut -d\  -f 3 | tr -d '\n' >> ../cputimes.dat 
		echo -n " " >> ../cputimes.dat
		grep 'matrix inverted' out_b0.log | tail -1 | cut -d\  -f 3 >> ../cputimes.dat
		cd ..
		let NPART=${NPART}+1
	done
	cd ..
	tri_mesh -r rstrt1_b0.grd refine.grd
	mv refine.grd rstrt1_b0.grd
	
	let NREFINE=${NREFINE}+1
done
	
cd ..

./make_plot.command

opendiff Baseline/ Results/
