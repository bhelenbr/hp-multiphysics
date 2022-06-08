#!/bin/bash

# This runs Mach 3 inviscid flow over the front half of a cylinder.
# 
# Downstream flow is initialized using incompressible potential flow over a cylinder at
# a velocity equal to that on the downstream side of the shock.
# 
# Make sure shock stabilization is turned on to get convergence.

# Generates a mesh for a curved bow shock
# where there are two domains before and after shock
#
# Uses cns which is based on primitive variables
# p, u, v, RT
# variables are non dimensional
# so p = 1/gamma
# u, v, = u/c
# RT = 1/gamma
# This makes rho = 1
# and c = 1

# Can't get this started with p=1 for some reason
# Have to start with p=4 and with the time_scheme: 4
# Doesn't start with time_scheme: 1

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

mpiexec -np 2 ${HP} run.inpt ${PETSC}

RESTART=$(mod_map -e run.inpt ntstep)
mod_map run.inpt ntstep 3
mod_map -c run.inpt auto_timestep_tries
mod_map run.inpt dtinv 0.0
mod_map run.inpt adapt 1
mod_map run.inpt error_estimator energy_norm
mod_map run.inpt error_target 1.0e-3
mod_map run.inpt restart ${RESTART}
mpiexec -np 2 ${HP} run.inpt ${PETSC}


#cd ..

#opendiff Results Baseline
