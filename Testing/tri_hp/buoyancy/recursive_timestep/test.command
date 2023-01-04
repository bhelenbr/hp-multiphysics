#!/bin/bash

# Runs a case with a facet point to test recursive timestepping
# Should not converge on time step 8 then do two finer steps
# last time step it fails

# Second case is to test the the case where 2 successful substeps
# are needed instead of 1 (nsuccesses 2)

# cd to the directory where the script resides.
# Necessary for platforms where script can be 
# double clicked from gui
cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/*/*/*/*}/bin
export PATH=${PATH}:${BINDIR}

# Make Results directory
if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm *

# copy input files into results directory
cp ../Inputs/* .

HP="mpiexec -np 2 ${MF} tri_hp_petsc run.inpt"

# generate mesh and remove unnecessary data files
tri_mesh generate
rm data*.grd

# These are the steps to run
# Get a reasonable steady temperature field
# Before starting the melting 
cp run.inpt startup.inpt
mod_map run.inpt b0_v1_hp_type plain
mod_map run.inpt b0_v2_hp_type plain
mod_map run.inpt b0_s3_hp_type inflow
mod_map run.inpt b0_s3_type symbolic
mod_map run.inpt b1_s3_hp_type dirichlet
mod_map run.inpt b1_s3_type symbolic
mod_map run.inpt b1_v1_hp_type plain
mod_map run.inpt b1_v2_hp_type plain
mod_map -c run.inpt mesh_movement
mod_map run.inpt ntstep 1
mod_map run.inpt dtinv 0.0
mod_map run.inpt restart_interval 1
# adapt has to be off so s3 doesn"t get messed up.
mod_map run.inpt adapt 0
${HP}
let RESTART=1

mv startup.inpt run.inpt
mod_map run.inpt ntstep 10
mod_map run.inpt recursive_timestep_levels 3
mod_map run.inpt restart $RESTART
mod_map run.inpt dtinv_prev 0.0
${HP}

mod_map run.inpt recursive_nsuccesses 2
${HP}
cd ..
opendiff Results/out_b0.log Baseline/out_b0.log
opendiff Results/out_b0.2.log Baseline/out_b0.2.log



