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
rm -rf *
cp ../Inputs/* .


let LOG2PMAX=3
let NGRIDMAX=6
let log2p=0

let ngrid=1
while [ $ngrid -lt ${NGRIDMAX} ]; do
	mod_map run.inpt b0_mesh square${ngrid}_b0.grd
	mod_map run.inpt ntstep 1
	let ngp=${ngrid}+1
	cp run.inpt square${ngrid}_b0.grd_bdry.inpt
	tri_mesh -r square${ngrid}_b0.grd square${ngp}_b0.grd
	let ngrid=${ngp}
done

mod_map run.inpt restart_interval 10000
while [ $log2p -lt ${LOG2PMAX} ]; do
	let ngrid=4-${log2p}
	let nmax=${NGRIDMAX}-${log2p}
	while [ ${ngrid} -le ${nmax} ]; do
		mkdir L2P_${log2p}_G${ngrid}
		cd L2P_${log2p}_G${ngrid}
		cp ../run.inpt .
		mod_map run.inpt ntstep $(echo "${ngrid}^4*2^(${log2p}*4)/8" | bc)
		mod_map run.inpt b0_mesh ../square${ngrid}_b0.grd
		mod_map run.inpt log2p ${log2p}
		${HP} run
		tail -2 output_b0.log | head -1 | cut -d\  -f2,4 >> ../cnvg${log2p}.dat
		grep L_2 output_b0.log | cut -d\  -f 2,4 > errs.dat
		cd ..
		let ngrid=${ngrid}+1
	done
	let log2p=${log2p}+1
done

cd ..

opendiff Results/ Baseline/
