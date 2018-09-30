#!/bin/bash
# Runs a case with a solidifying surface and periodic boundaries
# Interface is planar, analytic solution is known
# Tests different interface formulations
# one_sided does not work because non_sparse and match_jacobian
# at endpoints is not consistent with what is done along interface
# Beware that no cases work if log2p = 1.  Bizarre.
# also, this will not converge with an inflow b.c. because
# the pressure becomes unconstrained.

# To run_petsc test case comment 
#CNVG=true
NORMAL=true
#ONE_SIDED=true
SYMMETRIC=true
PRECONDITION=true

#PETSC="-stop_for_debugger"

cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/*/*/*/*}/bin
export PATH=${PATH}:${BINDIR}


if [ -e Results_petsc ]; then
	cd Results_petsc
else
	mkdir Results_petsc
	cd Results_petsc
fi
rm *

cp ../Inputs/* .

tri_mesh generate
rm data*.grd

HP="mpiexec -np 2 tri_hp_petsc run_petsc.inpt ${PETSC}"

if [ -n "$CNVG" ]; then
	let nrefinements=4
	let ngrid=0
	let ntstep=4
	while [ $ngrid -lt $nrefinements ]; do
		${HP} ${PETSC}
		grep L_2 out_b0.log | sed "s/\#L_2://g" | sed "s/L_inf \([-+.e0-9]*\)[ ]*[0-9]*/\1/g" > err.dat
		l2=$(awk 'BEGIN {max = 0} {if ($5>max) max=$5} END {print max}' err.dat )
		li=$(awk 'BEGIN {max = 0} {if ($6>max) max=$6} END {print max}' err.dat )
		echo "$l2 $li" >> cnvg.dat 

		tri_mesh -r rstrt1_b0.grd temp_b0.grd
		tri_mesh -r rstrt1_b1.grd temp_b1.grd
		mv temp_b0.grd rstrt1_b0.grd
		mv temp_b1.grd rstrt1_b1.grd
		let ntstep=${ntstep}+${ntstep}
		mod_map run_petsc.inpt ntstep ${ntstep}
		let ngrid=${ngrid}+1
	done
fi

if [ -n "$NORMAL" ]; then
	if [ -e NORMAL ]; then
		cd NORMAL
	else
		mkdir NORMAL
		cd NORMAL
	fi
	rm *

	cp ../*.grd .
	cp ../run_petsc.inpt .
	${HP}
	cd ..
fi

if [ -n "$ONE_SIDED" ]; then
	if [ -e ONE_SIDED ]; then
		cd ONE_SIDED
	else
		mkdir ONE_SIDED
		cd ONE_SIDED
	fi
	rm *

	cp ../*.grd .
	cp ../run_petsc.inpt .
	mod_map run_petsc.inpt b0_s3_one_sided 1
	${HP}
	if [ -n "$PRECONDITION" ]; then
		mod_map run_petsc.inpt b0_s3_precondition 1
		${HP}
	fi
	cd ..
fi

if [ -n "$SYMMETRIC" ]; then
	if [ -e SYMMETRIC ]; then
		cd SYMMETRIC
	else
		mkdir SYMMETRIC
		cd SYMMETRIC
	fi
	rm *

	cp ../*.grd .
	cp ../run_petsc.inpt .
	mod_map run_petsc.inpt b1_v1_hp_type hp_deformable_free_pnt
	mod_map run_petsc.inpt b1_v2_hp_type hp_deformable_free_pnt
	mod_map run_petsc.inpt b0_s3_symmetric 1
	${HP}
	if [ -n "$PRECONDITION" ]; then
		mod_map run_petsc.inpt b0_s3_precondition 1
		${HP}
	fi
	cd ..
fi
	
cd ..

opendiff Baseline_petsc/ Results_petsc/
