#!/bin/bash
# Runs a case that inverts the mass matrix
# using the approximate mass matrix
# Set #define DIRK 1 in tri_mesh/blocks.cpp

cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/*/*/*/*}/bin
export PATH=${PATH}:${BINDIR}

HP="tri_hp"


if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi

rm *
cp ../Inputs/* .

FULL_TEST=1

if [ -z ${FULL_TEST} ]; then
	let LOG2PMAX=3
	let NGRIDMAX=2
	let log2p=2
else
	let LOG2PMAX=3
	let NGRIDMAX=8
	let log2p=0
fi

while [ $log2p -lt ${LOG2PMAX} ]; do
	let ngrid=1
	let nedge=1
	while [ $ngrid -lt ${NGRIDMAX} ]; do
		mod_map run.inpt b0_mesh square${ngrid}_b0.grd
		mod_map run.inpt log2p ${log2p}
		mod_map run.inpt logfile approx${ngrid}.p${log2p}
		mod_map run.inpt ncycle 1
		${HP} run
		tail -2 approx${ngrid}.p${log2p}_b0.log | head -1 | cut -d\  -f2,4 >> cnvg${log2p}.dat
		
		mod_map run.inpt ncycle 40
		mod_map run.inpt logfile full${ngrid}.p${log2p}
		${HP} run
		tail -2 full${ngrid}.p${log2p}_b0.log | head -1 | cut -d\  -f2,4 >> cnvg${log2p}.dat
		
		let ngp=${ngrid}+1
		cp run.inpt square${ngrid}_b0.grd_bdry.inpt
		tri_mesh -r square${ngrid}_b0.grd square${ngp}_b0.grd
		let ngrid=${ngp}
		let nedge=2*${nedge}
	done
	let log2p=${log2p}+1
done

cd ..

./make_plot.command > Results/rates.dat

opendiff Results/ Baseline/
