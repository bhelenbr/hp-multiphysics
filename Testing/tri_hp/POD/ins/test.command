#!/bin/bash
# Runs oscillating flow over a sphere DNS
# Then generates modes and re-runs
# drag compared in plot

DNS=true
POD_GEN=true
POD_SIM=true
POD_SIM_PETSC=true

cd "$(dirname "$0")"

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi


if [ -n "$DNS" ]; then
	if [ -e DNS ]; then
		cd DNS
	else
		mkdir DNS
		cd DNS
	fi
	rm *
	
	cp ../../Inputs/* .

	mpiexec -np 1 tri_hp_axi_petsc run_steady.inpt
	rename rstrt1 rstrt0 rstrt1*

	mpiexec -np 1 tri_hp_axi_petsc run.inpt

	delete_data.bash 0 30

	# This is easy to understand but doesn't get the timestep numbers
	# grep -A1 pressure out_b0.log | grep "\[" | cut -d\  -f2-4 > drag.dat

	# This is a wee bit complicated but get the timestep number and the drag values
	grep -A1 'pressure\|TIMESTEP' out_b0.log | grep 'SUBSTEP: 2\|\[' | sed -e ':a' -e 'N' -e '$!ba' -e 's/2\n/ /g' | cut -d\  -f2,6-8 > drag.dat
	cd ..
fi


if [ -n "$POD_GEN" ]; then
	if [ -e pod_gen ]; then
		cd pod_gen
	else
		mkdir pod_gen
		cd pod_gen
	fi
	rm *

	cp ../DNS/rstrt* .
	cp ../../Inputs/* .

	tri_hp_axi POD_generate.inpt
	rm temp*
	rm rstrt*

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
	
	cp ../POD_GEN/mode* .
	cp ../DNS/rstrt45* .
	cp ../../Inputs/* .


	# Solve for the modes
	tri_hp_axi POD_simulate.inpt

	# This is easy to understand but doesn't get the timestep numbers
	# grep -A1 pressure out_b0.log | grep "\[" | cut -d\  -f2-4 > drag.dat

	# This is a wee bit complicated but get the timestep number and the drag values
	grep -A1 'pressure\|TIMESTEP' out_b0.log | grep ': 2\|\[' | sed -e ':a' -e 'N' -e '$!ba' -e 's/2\n/ /g' | cut -d\  -f2,6-8 > drag.dat

	cd ..
fi

if [ -n "$POD_SIM_PETSC" ]; then
	if [ -e pod_sim_petsc ]; then
		cd pod_sim_petsc
	else
		mkdir pod_sim_petsc
		cd pod_sim_petsc
	fi
	rm *
	
	cp ../POD_GEN/mode* .
	cp ../DNS/rstrt45* .
	cp ../../Inputs/* .


	# Solve for the modes
	mpiexec -np 1 tri_hp_axi_petsc POD_simulate.inpt

	# This is easy to understand but doesn't get the timestep numbers
	# grep -A1 pressure out_b0.log | grep "\[" | cut -d\  -f2-4 > drag.dat

	# This is a wee bit complicated but get the timestep number and the drag values
	grep -A1 'pressure\|TIMESTEP' out_b0.log | grep ': 2\|\[' | sed -e ':a' -e 'N' -e '$!ba' -e 's/2\n/ /g' | cut -d\  -f2,6-8 > drag.dat

	cd ..
	cp ../Inputs/make_plot.command .
	./make_plot.command
fi

cd ..

opendiff Baseline/ Results/
opendiff Results/pod_sim Results/pod_sim_petsc

