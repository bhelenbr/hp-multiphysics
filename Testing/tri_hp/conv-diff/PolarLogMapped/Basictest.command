#!/usr/bin/env bash
set -euo pipefail
#set -x
export FI_PROVIDER=tcp

HP="mpiexec -np 1 tri_hp_petsc"
# Testing accuracy for a case with a singular point

epsVal="$1"   # required argument

# cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/Testing/*}/bin
export PATH=${BINDIR}:${PATH}

# mkdir -p Results
# cd Results
# rm -rf ./*

# cp ../Inputs/* .

mod_map generate.inpt eps "$epsVal"
tri_mesh generate.inpt

cp generate.inpt run.inpt

mod_map run.inpt b0_mesh rstrt1_b0.grd
mod_map run.inpt logfile output
mod_map run.inpt adapt 0
mod_map run.inpt ncycle 10
mod_map run.inpt ntstep 1

extract_dof() {
    awk '/DOF:/ {print $6}' output_b0.log
}

append_metrics() {
    printf "%s" "$(extract_dof)" >> cnvg.dat
    tail -3 output_b0.log | awk 'NR==1 {print " "$2" "$4}' >> cnvg.dat
}


log2p=0

while (( log2p < 3 )); do
    
    workdir="log2p${log2p}"
    mkdir "$workdir"
    cp run.inpt rstrt1_b0.grd "$workdir/"
    cd "$workdir"

    mod_map run.inpt log2p "$log2p"
	$HP run.inpt

	append_metrics
	
	ngrids=4
	ngrid=1
	restart=1

	while (( ngrid < ngrids )); do
		# Refine solution
		mod_map run.inpt refineby2 1
		mod_map run.inpt adapt 1
		mod_map run.inpt restart "$restart"
		$HP run.inpt
		((restart++))
		
		# Run case using restart file
		mod_map run.inpt restart "$restart"
		mod_map run.inpt adapt 0
		mod_map run.inpt refineby2 0
		mod_map run.inpt b0_mesh "rstrt${restart}_b0.nc"
		$HP run.inpt
		((restart++))
		
        append_metrics
        ((ngrid++))
	done
	cd ..
	((log2p++))
done
cd ..
#./make_plot.command > Results/rates.dat
#opendiff Results/ Baseline/