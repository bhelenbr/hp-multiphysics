#!/bin/bash
# Rigid body translation using the moving-mesh formulation
# No gravity translating periodic flow
# There are several ways to set-up the communication
# Currently the only one that works is phased communication
# Solution should be a uniform flow and a constant velocity
# translation of sinusoidal interface

cd "$(dirname "$0")"

# To run_mgrid test case comment 
#CNVG=true
#MGRID_PHASED_WALLS_FIRST=true
#MGRID_PHASED_SURFACE_FIRST=true
#MGRID_ONEPASS=true

cd "$(dirname "$0")"

HP="tri_hp"
MESH="tri_mesh"


if [ -e Results_mgrid ]; then
	cd Results_mgrid
else
	mkdir Results_mgrid
	cd Results_mgrid
fi
rm -rf *

cp ../Inputs/* .

cp run_mgrid.inpt generate.inpt
mod_map generate.inpt logfile generate
mod_map generate.inpt restart_interval 1
mod_map generate.inpt b0_filetype 8
mod_map generate.inpt b1_filetype 8
mod_map generate.inpt b0_growth\ factor 1000.0
mod_map generate.inpt b1_growth\ factor 1000.0
mod_map generate.inpt b0_s5_type symbolic_comm
mod_map generate.inpt b1_s5_type symbolic_comm
mod_map generate.inpt ngrid 1
mod_map generate.inpt log2p 0
mod_map generate.inpt extra_finest_levels 0
mod_map generate.inpt adapt 1
mod_map generate.inpt ntstep 1
mod_map -c generate.inpt ncycle

if [ -n "$CNVG" ]; then
	NGRID=5
	mod_map generate.inpt RES 1.0
	tri_mesh generate
	mv rstrt1_b0.grd bot1.grd
	mv rstrt1_b1.grd top1.grd
	rm data*
	mod_map run_mgrid.inpt b0_s5_report 1

	let ngrid=1
	let dtinv=4
	while [ $ngrid -lt ${NGRID} ]; do
		let dtinv=${dtinv}*2
		let ntstep=4*dtinv
		mkdir grid${ngrid}
		cd grid${ngrid}
		cp ../run_mgrid.inpt .
		cp ../*.grd .
		mod_map run_mgrid.inpt ngrid ${ngrid}
		mod_map run_mgrid.inpt dtinv ${dtinv}
		if [ -n "$CNVG" ]; then
			mod_map run_mgrid.inpt ntstep ${ntstep}
		fi
		${HP} run_mgrid
		grep b0_s5\ l2 out_b0.log | cut -d\    -f4 > l2_iface${ngrid}.dat
		mv l2_iface${ngrid}.dat ..
		let ngp=${ngrid}+1
		cd ..
		${MESH} -r bot1 bot1
		${MESH} -r top1 top1
		let ngrid=${ngp}
	done
	cd ..
else
	mod_map generate.inpt RES 0.25
	tri_mesh generate
	mv rstrt1_b0.grd bot1.grd
	mv rstrt1_b1.grd top1.grd
	rm data*
	${HP} run_mgrid
	cd ..
	opendiff Baseline_mgrid Results_mgrid
fi

