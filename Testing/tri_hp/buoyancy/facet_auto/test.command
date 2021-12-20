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
# -stop_for_debugger"

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
if [ "$?" -ne "0" ]; then
  echo "exited from initial T solution"
  exit 1
fi 
let RESTART=1

mv startup.inpt run.inpt
# 100 steps with small time step
mod_map run.inpt auto_timestep 1
mod_map run.inpt auto_timestep_ratio 1.25
mod_map run.inpt ntstep 100
mod_map run.inpt restart $RESTART
mod_map run.inpt restart_interval 1
mod_map run.inpt dtinv_prev 0.0
#mod_map run.inpt debug_output 1
#mod_map run.inpt under_relaxation 0.66
${HP}
if [ "$?" -ne "0" ]; then
  echo "exited from initial T solution"
  exit 1
fi 
RESTART=100
echo "unsteady evolution 1 ${RESTART}"

# Steady State
mod_map -c run.inpt auto_timestep
mod_map run.inpt restart ${RESTART}
mod_map run.inpt ntstep 1
mod_map run.inpt restart_interval 1
mod_map -u run.inpt under_relaxation
mod_map -d run.inpt dtinv_prev
mod_map run.inpt dtinv 0.0
${HP}
if [ "$?" -ne "0" ]; then
  echo "exited from initial steady-state ${RESTART}"
  exit 1
fi 

exit 1

let RESTART=${RESTART}+1
echo "initial steady-state ${RESTART}"
# 501

# MAKE K2DN_MAX infinite
mod_map run.inpt restart ${RESTART}
mod_map run.inpt b0_s3_K2Dn_max 1e10
${HP}
if [ "$?" -ne "0" ]; then
	echo "exited from changing K2DN_max"
	exit 1
fi 

let RESTART=${RESTART}+1
echo "changing K2DN_max ${RESTART}"
# 502

# Get high order solution (try to get 1step steady-state first)
mod_map run.inpt dtinv 0.0
mod_map run.inpt log2p 1
mod_map run.inpt restart ${RESTART}
mod_map run.inpt adapt 0
mod_map run.inpt under_relaxation 0.25
${HP}
if [ "$?" -ne "0" ]; then
	mod_map run.inpt dtinv "1e4*(dtinv1 +dtinv3)"
	mod_map run.inpt ntstep 20
	mod_map -c run.inpt under_relaxation
	${HP}
	if [ "$?" -ne "0" ]; then
		echo "exited from log2p = 1"
		exit 1
	fi 
	
	let RESTART=${RESTART}+20
	mod_map run.inpt dtinv 0.0
	mod_map run.inpt restart ${RESTART}
	mod_map run.inpt ntstep 1
	mod_map -u run.inpt under_relaxation
	${HP}
	if [ "$?" -ne "0" ]; then
		echo "exited from log2p = 1"
		exit 1
	fi
fi 	

let RESTART=${RESTART}+1
echo "log2p = 1  ${RESTART}"

# Get high order solution (try to get 1step steady-state first)
mod_map run.inpt log2p 2
mod_map run.inpt restart ${RESTART}
mod_map run.inpt adapt 0
mod_map run.inpt dtinv 0.0
mod_map run.inpt ntstep 1
mod_map -u run.inpt under_relaxation
${HP}
if [ "$?" -ne "0" ]; then
	mod_map run.inpt dtinv "1e4*(dtinv1 +dtinv3)"
	mod_map run.inpt ntstep 20
	mod_map -c run.inpt under_relaxation
	${HP}
	if [ "$?" -ne "0" ]; then
		echo "exited from log2p = 2"
		exit 1
	fi
	
	let RESTART=${RESTART}+20
	mod_map run.inpt dtinv 0.0
	mod_map run.inpt restart ${RESTART}
	mod_map run.inpt ntstep 1
	mod_map -u run.inpt under_relaxation
	${HP}
	if [ "$?" -ne "0" ]; then
		echo "exited from log2p = 1"
		exit 1
	fi
fi
	

let RESTART=${RESTART}+1
echo "log2p = 2 ${RESTART}"
mod_map run.inpt under_relaxation 1.0

# Get steady-state
mod_map run.inpt dtinv 0.0
mod_map run.inpt restart ${RESTART}
mod_map run.inpt ntstep 1
${HP}
if [ "$?" -ne "0" ]; then
	echo "exited from log2p = 2 steady state"
	exit 1
fi

let RESTART=${RESTART}+1
echo "log2p = 2 steady state ${RESTART}"

# Turn on Mesh adaptation
mod_map -u run.inpt error_estimator
mod_map run.inpt adapt 1
mod_map run.inpt restart ${RESTART}
mod_map run.inpt ntstep 3
${HP}
if [ "$?" -ne "0" ]; then
	echo "exited from initial mesh adaptation"
	exit 1
fi

let RESTART=${RESTART}+3
echo "mesh adaptation ${RESTART}"

################################################
# SX IS THE PULL SPEED
# NOW SLOWLY INCREASE SX
################################################
mod_map -c run.inpt under_relaxation
mod_map run.inpt dtinv 0.0
mod_map run.inpt ntstep 1
delete_series data 0 1
delete_series rstrt 2 1

SX="$(mod_map -e run.inpt sx | cut -d. -f1)0"
echo ${SX} > sx.dat
head -2 "data${RESTART}_b1.dat" | tail -1 | cut -d\  -f1 > xle.dat

SXMAX=1000000
FACTOR="1.05"
DSX=$(echo "(${FACTOR}-1)*${SX}"| bc -l | cut -d. -f1)
let ATTEMPTS=0
let MAXATTEMPTS=3

while [ ${ATTEMPTS} -lt ${MAXATTEMPTS} ] && [ $(echo "${SX} < ${SXMAX}" | bc -l) -eq 1 ]; do
	SXOLD=${SX}
	let SX=${SX}+${DSX}
	echo "SX ${SX}"
	
	mod_map run.inpt sx "${SX}e-7/(d0*tsi)"
	mod_map run.inpt restart ${RESTART}

	${HP}
	if [ "$?" -ne "0" ]; then
		rm core*
		rm neg*
		rm abort*
		SX=${SXOLD}
		DSX=$(echo "${DSX}/4"| bc -l | cut -d. -f1)
		let ATTEMPTS=${ATTEMPTS}+1
	  	mod_map run.inpt extrapolate 0.25
	else
		# TEST TO MAKE SURE CONVERGED
		RES=$(grep ^[1-9] out_b0.log | tail -1 | cut -d\  -f 3)
		RES=${RES/[eE]/*10^\(}\)
		RES=${RES/\([+]/\(}
		CNVG=$(echo "${RES} < 1.0*10^(-5)" | bc -l)
		if [ "$CNVG" != "1" ]; then
			DSX=$(echo "${DSX}/4"| bc -l | cut -d. -f1)
			let ATTEMPTS=${ATTEMPTS}+1
	  		mod_map run.inpt extrapolate 0.25
		else
			let RESTART=${RESTART}+1
			mod_map run.inpt extrapolate 1.0
			echo ${SX} >> sx.dat
			head -2 "data${RESTART}_b1.dat" | tail -1 | cut -d\  -f1 >> xle.dat
		fi
	fi
done
../plot.py

rm core* abort* net* rstrt*.nc
