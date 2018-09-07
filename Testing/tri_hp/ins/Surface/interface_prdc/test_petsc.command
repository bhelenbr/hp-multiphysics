#!/bin/bash
# Rigid body translation using the moving-mesh formulation
# No gravity translating periodic flow
# There are several ways to set-up the communication
# Currently the only one that works is phased communication
# Solution should be a uniform flow and a constant velocity
# translation of sinusoidal interface

# phased communication does not really work with petsc
# for this to work, J_mpi would need to be message passed
# nnzero_mpi would also need to be message passed
# this case works only if you do message passing to
# local block first (periodic boundaries)

cd "$(dirname "$0")"

#PETSC_FLAGS="-stop_for_debugger"

# To run_petsc test case comment 
#CNVG=true
PETSC_PHASED_WALLS_FIRST=true
#PETSC_PHASED_SURFACE_FIRST=true
PETSC_ONEPASS=true

cd "$(dirname "$0")"

HP="mpiexec -np 2 tri_hp_petsc"
MESH="tri_mesh"


if [ -e Results_petsc ]; then
	cd Results_petsc
else
	mkdir Results_petsc
	cd Results_petsc
fi
rm *

cp ../Inputs/* .

cp run_petsc.inpt generate.inpt
mod_map generate.inpt nblock 2
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
	mod_map run_petsc.inpt b0_s5_report 1

	let ngrid=1
	let dtinv=4
	while [ $ngrid -lt ${NGRID} ]; do
		let dtinv=${dtinv}*2
		let ntstep=4*dtinv
		mkdir grid${ngrid}
		cd grid${ngrid}
		cp ../run_petsc.inpt .
		cp ../*.grd .
		mod_map run_petsc.inpt ngrid ${ngrid}
		mod_map run_petsc.inpt dtinv ${dtinv}
		if [ -n "$CNVG" ]; then
			mod_map run_petsc.inpt ntstep ${ntstep}
		fi
		${HP} run_petsc
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
	
	if [ -n "$PETSC_PHASED_WALLS_FIRST" ]; then
		if [ -e PETSC_PHASED_WALLS_FIRST ]; then
			cd PETSC_PHASED_WALLS_FIRST
		else
			mkdir PETSC_PHASED_WALLS_FIRST
			cd PETSC_PHASED_WALLS_FIRST
		fi
		rm *
	
		cp ../*.grd .
		cp ../run_petsc.inpt .
		${HP} run_petsc ${PETSC_FLAGS}
		cd ..
	fi
	
	if [ -n "$PETSC_PHASED_SURFACE_FIRST" ]; then
		if [ -e PETSC_PHASED_SURFACE_FIRST ]; then
			cd PETSC_PHASED_SURFACE_FIRST
		else
			mkdir PETSC_PHASED_SURFACE_FIRST
			cd PETSC_PHASED_SURFACE_FIRST
		fi
		rm *
	
		cp ../*.grd .
		cp ../run_petsc.inpt .
		mod_map -c run_petsc.inpt b0_s5_phase1
		mod_map -c run_petsc.inpt b1_s5_phase1
		mod_map -u run_petsc.inpt b0_s2_phase1
		mod_map -u run_petsc.inpt b1_s1_phase1
		${HP} run_petsc ${PETSC_FLAGS}
		cd ..
	fi
	
	if [ -n "$PETSC_ONEPASS" ]; then
		if [ -e PETSC_ONEPASS ]; then
			cd PETSC_ONEPASS
		else
			mkdir PETSC_ONEPASS
			cd PETSC_ONEPASS
		fi
		rm *
	
		cp ../*.grd .
		cp ../run_petsc.inpt .
		mod_map -c run_petsc.inpt b0_s5_phase1
		mod_map -c run_petsc.inpt b1_s5_phase1
		mod_map -u run_petsc.inpt b0_v1_type
		mod_map -u run_petsc.inpt b1_v1_type
		${HP} run_petsc ${PETSC_FLAGS}
		cd ..
	fi
	
	cd ..
	opendiff Baseline_petsc Results_petsc
fi

