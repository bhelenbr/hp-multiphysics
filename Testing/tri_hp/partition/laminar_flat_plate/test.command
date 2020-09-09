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

tri_mesh generate.inpt

# HP="mpiexec -np 1 tri_hp_petsc -stop_for_debugger"

#PETSC="-stop_for_debugger"

mpiexec -np 1 ${HP} run.inpt ${PETCSC}
echo -n "1 " >> cputimes.dat
tail -1 out_b0.log | cut -d \  -f 3 >> cputimes.dat
tail -1 out_b0.log | cut -d \  -f 3 | tr -d '\n' >> cputimes.dat
echo -n " " >> cputimes.dat
grep 'jacobian made' out_b0.log | tail -1 | cut -d\  -f 3 | tr -d '\n' >> cputimes.dat 
echo -n " " >> cputimes.dat
grep 'matrix inverted' out_b0.log | tail -1 | cut -d\  -f 3 >> cputimes.dat

let NPART=2
let NPROC=4

while [ ${NPART} -le ${NPROC} ]; do
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
	mpiexec -np ${NPART} ${HP} partition.inpt ${PETSC}
	echo -n "${NPART} "  >> ../cputimes.dat
	tail -1 out_b0.log | cut -d \  -f 3 | tr -d '\n' >> ../cputimes.dat
	echo -n " " >> ../cputimes.dat
	grep 'jacobian made' out_b0.log | tail -1 | cut -d\  -f 3 | tr -d '\n' >> ../cputimes.dat 
	echo -n " " >> ../cputimes.dat
	grep 'matrix inverted' out_b0.log | tail -1 | cut -d\  -f 3 >> ../cputimes.dat
	cd ..
	let NPART=${NPART}+1
done
../plot.py

cd ..

opendiff Baseline/ Results/
