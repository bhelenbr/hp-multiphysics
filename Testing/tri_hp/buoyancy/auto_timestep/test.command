#!/bin/bash

# Runs a case with a facet point
# Flow is uniform flow
# All cases except symmetric perform similarly
# under_relaxation is needed to get log2p=0 to log2p=1 switch to work
# once it worked without it, but I never figured out what was different
# uncomment full_test to do the whole thing
# otherwise does one time step only
# For the short test if adapt is on, it will diverge at log2p = 2
# with adapt off it will converge (with no under_relaxation)
# symmetric has not made it all the way to the turning point
# Because v,p = 0, can get poor convergence because of numerical Jacobian
# Sensitive to eps_a in evaluating Jacobian when using dw = dw*eps_r +eps_a
# With new kinetic expression it stopped one time step before when using the old kinetic expression!

# This cases is highly sensitive to the LU factorization.  I got the best results using
# mumps, with petsc version 3.7.7.  After that, it got worse

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

FULL_TEST=1

HP="mpiexec -np 2 ${MF} tri_hp_petsc run.inpt"

# generate mesh and remove unnecessary data files
tri_mesh generate
rm data*.grd

# Various ways of running
#mod_map run.inpt b0_s3_one_sided 1
#mod_map run.inpt b0_s3_precondition 1

# This way never made it from log2p = 1 to log2p = 2
# mod_map run.inpt b0_s3_symmetric 1
# mod_map run.inpt b1_v1_hp_type hp_deformable_free_pnt
# mod_map run.inpt b1_v2_hp_type melt_facet_pt
	
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
# 6 steps with small time step 1/2 time step after step 4
mod_map run.inpt log2p 1
mod_map run.inpt ntstep 4
mod_map run.inpt auto_timestep 1
mod_map run.inpt restart $RESTART
mod_map run.inpt restart_interval 3
mod_map run.inpt dtinv_prev 0.0
#mod_map run.inpt debug_output 1
#mod_map run.inpt under_relaxation 0.66
${HP}

mod_map run.inpt restart 3
mod_map -c run.inpt auto_timestep
mod_map run.inpt dtinv_prev "1e4*(dtinv1 +dtinv3)"
mod_map run.inpt dtinv "0.5e4*(dtinv1+dtinv3)"
mod_map run.inpt ntstep 1
${HP}

