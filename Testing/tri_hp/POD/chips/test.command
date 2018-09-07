#!/bin/bash
# Runs a cd POD case where there are 2 POD blocks
# separated by a pod boundary that is composed of
# two separate boundaries.

SOLVE=true
POD_GEN=true
POD_SIM=true

HP=tri_hp

cd "$(dirname "$0")"

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi

if [ -n "$SOLVE" ]; then
	if [ -e solve ]; then
		cd solve
	else
		mkdir solve
		cd solve
	fi
	rm *

	# Solve for the modes
	cp ../../Inputs/blocks*.inpt .
	cp ../../Inputs/*.d .
	tri_mesh blocks_mesh
	rm data*

	# adapt mesh
	${HP} blocks_adapt
	
	# perform simulations
	rename rstrt5 initial rstrt5_*.grd
	rm rstrt*
	${HP} blocks_solve

	cd ../
fi

if [ -n "$POD_GEN" ]; then
	if [ -e pod_gen ]; then
		cd pod_gen
	else
		mkdir pod_gen
		cd pod_gen
	fi
	rm *

	cp ../solve/rstrt* .
	cp ../solve/initial* .

	# Solve for the modes
	cp ../../Inputs/pod_gen* .
	#${HP} pod_gen1block
	${HP} pod_gen
	cd ..
fi

if [ -n "$POD_SIM" ]; then
	if [ -e pod_sim ]; then
		cd pod_sim
	else
		mkdir pod_sim
		cd pod_sim
	fi
	rm *
	
	cp ../pod_gen/mode* .
	cp ../pod_gen/initial* .
	cp ../pod_gen/coeff1* .

	# Solve for the modes
	cp ../../Inputs/pod_sim* .
	#${HP} pod_sim1block
	${HP} pod_sim
	cd ..
fi

cd ..

opendiff Baseline/ Results/
