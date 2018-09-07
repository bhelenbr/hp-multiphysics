#!/bin/bash
# Rigid body translation using the level-set formulation
# # To tun this case, must make sure flow is not updated (u=1)
# Have to brute force make sure in level-set tstep.cpp: setup_preconditioner
# flow is not updated
# If you want reinitilization have to make all reinit_bc's inactive
# then turn on communication in reinit_minvrt

cd "$(dirname "$0")"

# To run test case comment 
CNVG=true

if [ -n "$CNVG" ]; then
	NGRID=4
else
	NGRID=2
fi


cd "$(dirname "$0")"

if [ -e $HOME/Codes/tri_hp/build/Release/tri_hp ]; then
  HP="$HOME/Codes/tri_hp/build/Release/tri_hp"
else
  HP="$HOME/bin/tri_hp"
fi
MESH="$HOME/bin/tri_mesh"


if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm -rf *

cp ../Inputs/* .

# GENERATE NESTED GRID REFINEMENT
cp run.inpt generate.inpt
mod_map generate.inpt ntstep 2
mod_map generate.inpt adapt 1
mod_map generate.inpt b0_filetype 8
mod_map generate.inpt b0_growth\ factor 1000.0
mod_map generate.inpt ngrid 1
mod_map generate.inpt log2p 0
mod_map generate.inpt extra_finest_levels 0
mod_map generate.inpt b0_mesh ./conv
${MESH} generate
cp data2_b0.grd conv.grd

let ngrid=1
let dtinv=4
while [ $ngrid -le ${NGRID} ]; do
	let dtinv=${dtinv}*2
	let ntstep=4*dtinv
	mkdir grid${ngrid}
	cd grid${ngrid}
	cp ../run.inpt .
	cp ../*.grd .
	mod_map run.inpt ngrid ${ngrid}
	mod_map run.inpt dtinv ${dtinv}
	if [ -n "$CNVG" ]; then
		mod_map run.inpt ntstep ${ntstep}
	fi
	${HP} run
	grep Contour out_b0.log | cut -d\    -f3 > l2_zero${ngrid}.dat
	mv l2_zero${ngrid}.dat ..
	let ngp=${ngrid}+1
	cd ..
	${MESH} -r conv conv
	let ngrid=${ngp}
done

cd ..