#!/bin/bash
# Runs a startup simulation for flow over a NACA0012
# Tests POD with partitioned simulations
# Also tests using the petsc Jacobian instead of 
# numerically evaluating the Jacobian.  This DOES NOT WORK!
# The communication in matchjacobian_rcv messes it up.
# to get this to work you have to turn off the diagonal swapping
# and the reduction of the jacobian row by a factor of 1/nmatch.
# not sure what I am going to do about that.

DNS=true
POD_GEN=true
POD_SIM=true
POD_SIM_PETSC=true

HP="tri_hp_petsc"

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
	
	tri_mesh generate.inpt
	
	mod_map run.inpt nblock 4
	mod_map run.inpt partition 1
	mod_map run.inpt logfile partition
	tri_partition run.inpt	
	
	mod_map partition.inpt logfile output
	mpiexec -np 4 ${HP} partition.inpt

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
	cp ../DNS/partition* .
	mod_map partition.inpt snapshots 5
	mod_map partition.inpt blocktype pod_ins_gen
	mpiexec -np 4 ${HP} partition.inpt
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
	cp ../DNS/partition* .
	cp ../DNS/rstrt1_* .
	mod_map partition.inpt restart 1
	mod_map partition.inpt nmodes 5
	mod_map partition.inpt ntstep 4
	mod_map partition.inpt blocktype pod_ins_sim

	# Simulate
	mpiexec -np 4 tri_hp_mpi partition.inpt
	
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
	cp ../DNS/partition* .
	cp ../DNS/rstrt1_* .
	mod_map partition.inpt restart 1
	mod_map partition.inpt nmodes 5
	mod_map partition.inpt ntstep 4
	mod_map partition.inpt blocktype pod_ins_sim
	mod_map partition.inpt petsc "-ksp_type preonly"

	# Simulate
	mpiexec -np 4 tri_hp_petsc partition.inpt
	
	cd ..
fi



cd ..

opendiff Results/pod_sim Results/pod_sim_petsc
opendiff Baseline/ Results/

