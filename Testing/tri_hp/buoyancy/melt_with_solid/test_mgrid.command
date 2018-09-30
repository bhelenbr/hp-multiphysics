#!/bin/bash
# Runs a case with a solidifying surface
# Interface is planar, analytic solution is known
# This only works for log2p=0.  I put in stuff for side modes,
# but it is not working.  Not sure how to combine the flow equations
# with the interface equations with higher order modes

# To run_petsc test case comment 
#CNVG=true
NORMAL=true

cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/*/*/*/*}/bin
export PATH=${PATH}:${BINDIR}

if [ -e Results_mgrid ]; then
	cd Results_mgrid
else
	mkdir Results_mgrid
	cd Results_mgrid
fi
rm *

cp ../Inputs/* .

tri_mesh generate
rm data*.grd

#HP="mpiexec -np 2 tri_hp_mpi -stop_for_debugger run_mgrid.inpt"
HP="mpiexec -np 2 tri_hp_mpi run_mgrid.inpt"


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
		mod_map run_mgrid.inpt ntstep ${ntstep}
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
	cp ../run_mgrid.inpt .
	${HP}
	cd ..
fi
	
cd ..

opendiff Baseline_mgrid/ Results_mgrid/
