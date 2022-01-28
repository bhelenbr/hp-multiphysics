#!/bin/bash

# This runs Mach 10 viscous flow over a full cylinder.
#
# Generates a mesh for a curved bow shock
# where there are two domains before and after shock
#
# Uses cns which is based on primitive variables
# p, u, v, RT
# variables are non dimensional
# so p = gamma
# u, v, = u/c
# RT = 1/gamma
# This makes rho = 1
# and c = 1
#
# Process:
# 
# 1) Run run.inpt until steady state is reached (log2p = 0, shock is stationary).
# 2) run run.inpt with log2p1 until steady state is reached (shock is now moving).
# 3) run run.inpt with log2p2 (shock is still moving).
# 4) Do mesh refinement.
# 
# Make sure shock stabilization is turned on for convergence.

cd "$(dirname "$0")"
BINDIR=${PWD%/*/*/*/*}/bin
export PATH=${PATH}:${BINDIR}

HP="tri_hp_petsc"
#PETSC="-stop_for_debugger"

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm *

cp ../Inputs/* .

tri_mesh generate

# To run flow only
cp run.inpt startup.inpt
mod_map startup.inpt b0_s3_type symbolic
mod_map startup.inpt b0_s3_hp_type outflow_supersonic
mod_map startup.inpt b1_s3_hp_type inflow
mod_map startup.inpt b1_s3_type symbolic
mod_map -c startup.inpt b0_v1_hp_type
mod_map -c startup.inpt b0_v2_hp_type
mod_map -c startup.inpt b1_v1_hp_type
mod_map -c startup.inpt b1_v2_hp_type
mpiexec -np 2 ${HP} startup.inpt ${PETSC}
RESTART=$(mod_map -e run.inpt ntstep)

mod_map run.inpt restart ${RESTART}
mpiexec -np 2 ${HP} run.inpt ${PETSC}
let RESTART=${RESTART}+$(mod_map -e run.inpt ntstep)

mod_map run.inpt restart ${RESTART}
mod_map run.inpt log2p 1
mpiexec -np 2 ${HP} run.inpt ${PETSC}
let RESTART=${RESTART}+$(mod_map -e run.inpt ntstep)

mod_map run.inpt restart ${RESTART}
mod_map run.inpt log2p 2
mpiexec -np 2 ${HP} run.inpt ${PETSC}

#cd ..

#opendiff Results Baseline
