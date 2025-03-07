#!/bin/bash
# Runs a case with a facet, a liquid free-surface, a solid translating surface
# These meet at a triple point that moves according to the growth angle

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

cp run.inpt generate.inpt
mod_map generate.inpt b0_mesh geometry_b0.d
mod_map generate.inpt b1_mesh geometry_b1.d
mod_map generate.inpt "growth factor" 4000
mod_map generate.inpt b0_s3_type symbolic_comm
mod_map generate.inpt b1_s3_type symbolic_comm
mod_map generate.inpt b1_s7_type symbolic
mod_map generate.inpt b0_s4_type symbolic
mod_map generate.inpt nblock 2
mod_map generate.inpt restart_interval 1
mod_map generate.inpt logfile generate
mod_map generate.inpt ncycle 0

# generate mesh and remove unnecessary data files
tri_mesh generate
rm data*.grd

# These are the steps to run
# Get a reasonable steady temperature field
# Before starting the melting 
# mod_map -d  means delete that line from the input file
# mod_map run.inpt keyword value ## Changes the value in the input file
cp run.inpt startup.inpt

mod_map run.inpt b0_v1_hp_type plain
mod_map run.inpt b0_v2_hp_type hp_deformable_fixed_pnt
mod_map run.inpt b1_v1_hp_type plain
mod_map run.inpt b1_v2_hp_type plain
mod_map run.inpt b1_v4_hp_type plain
mod_map run.inpt b0_s3_hp_type inflow
mod_map run.inpt b0_s3_type symbolic
mod_map run.inpt b1_s3_type symbolic
mod_map run.inpt b1_s3_hp_type dirichlet
mod_map run.inpt b1_s7_type symbolic
mod_map run.inpt b1_s7_hp_type plain
mod_map run.inpt b1_s7_flux0 convflux+radiation*epss

mod_map run.inpt 
mod_map run.inpt ntstep 1
mod_map run.inpt dtinv 0.0
mod_map run.inpt restart_interval 1
${HP}
if [ "$?" -ne "0" ]; then
  echo "exited from initial T solution"
  exit 1
fi 
let RESTART=0
rename rstrt1 rstrt0 rstrt*.nc
rm data*

if [ -z ${FULL_TEST} ]; then
	mv startup.inpt run.inpt
	mod_map run.inpt ntstep 1
	mod_map run.inpt restart ${RESTART}
	${HP}
	if [ "$?" -ne "0" ]; then
	  echo "exited from unsteady evolution 1"
	  exit 1
	fi
	rm rstrt*.bin
	cd .. 
	opendiff Baseline/ Results/
	exit 0
fi

mv startup.inpt run.inpt
# 100 steps with small time step
mod_map run.inpt ntstep 100
mod_map run.inpt restart ${RESTART}
mod_map run.inpt restart_interval 100
${HP}
if [ "$?" -ne "0" ]; then
  echo "exited from unsteady evolution 1"
  exit 1
fi 
let RESTART=100
echo "unsteady evolution 1 $RESTART"

# More steps with larger time step
mod_map run.inpt restart ${RESTART}
mod_map run.inpt dtinv_prev "1e4*(dtinv1 +dtinv3)"
mod_map run.inpt dtinv "1e3*(dtinv1+dtinv3)"
${HP}
if [ "$?" -ne "0" ]; then
  echo "exited from unsteady evolution 3"
  exit 1
fi
let RESTART=${RESTART}+100
echo "unsteady evolution 3 $RESTART"
# 300


# More steps with larger time step
mod_map run.inpt restart ${RESTART}
mod_map run.inpt dtinv_prev "1e3*(dtinv1 +dtinv3)"
mod_map run.inpt dtinv "1e2*(dtinv1+dtinv3)"
${HP}
if [ "$?" -ne "0" ]; then
  echo "exited from unsteady evolution 4"
  exit 1
fi
let RESTART=${RESTART}+100
echo "unsteady evolution 4 $RESTART"
# 400

# More steps with larger time step
mod_map run.inpt restart ${RESTART}
mod_map run.inpt dtinv_prev "1e2*(dtinv1 +dtinv3)"
mod_map run.inpt dtinv "1e1*(dtinv1+dtinv3)"
${HP}
if [ "$?" -ne "0" ]; then
  echo "exited from unsteady evolution 5"
  exit 1
fi
let RESTART=${RESTART}+100
echo "unsteady evolution 5 $RESTART"
# 500

# Steady State
mod_map run.inpt restart ${RESTART}
mod_map run.inpt ntstep 1
mod_map run.inpt restart_interval 1
mod_map -u run.inpt under_relaxation
mod_map -d run.inpt dtinv_prev
mod_map run.inpt dtinv "0.0"
${HP}
if [ "$?" -ne "0" ]; then
  echo "exited from initial steady-state"
  exit 1
fi
let RESTART=${RESTART}+1
echo "initial steady-state $RESTART"
# 401

# MAKE K2DN_MAX infinite
mod_map run.inpt restart ${RESTART}
mod_map run.inpt b0_s3_K2Dn_max 1e10
${HP}
if [ "$?" -ne "0" ]; then
  echo "exited from changing K2DN_max"
  exit 1
fi
let RESTART=${RESTART}+1
echo "changing K2DN_max $RESTART"
# 401

# Get high order solution (try to get 1step steady-state first)
mod_map run.inpt dtinv "0.0"
mod_map run.inpt log2p 1
mod_map run.inpt restart ${RESTART}
mod_map run.inpt adapt 0
mod_map run.inpt under_relaxation 0.25
${HP}
if [ "$?" -ne "0" ]; then
	mod_map run.inpt dtinv "1e4*(dtinv1 +dtinv3)"
	mod_map run.inpt ntstep 100
	mod_map -c run.inpt under_relaxation
	${HP}
	if [ "$?" -ne "0" ]; then
	  echo "exited from log2p = 1"
	  exit 1
	fi
	let RESTART=${RESTART}+20
	mod_map run.inpt dtinv "0.0"
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
echo "log2p = 1 ${RESTART}"

# Get high order solution (try to get 1step steady-state first)
mod_map run.inpt log2p 2
mod_map run.inpt restart ${RESTART}
mod_map run.inpt adapt 0
mod_map run.inpt dtinv "0.0"
mod_map run.inpt ntstep 1
mod_map -u run.inpt under_relaxation
${HP}
if [ "$?" -ne "0" ]; then
	mod_map run.inpt dtinv "1e4*(dtinv1 +dtinv3)"
	mod_map run.inpt ntstep 100
	mod_map -c run.inpt under_relaxation
	${HP}
	if [ "$?" -ne "0" ]; then
	  echo "exited from log2p = 2"
	  exit 1
	fi
	let RESTART=${RESTART}+20
	mod_map run.inpt dtinv "0.0"
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
delete_series data 0 1 
delete_series rstrt 2 1


################################################
# SX IS THE PULL SPEED
# NOW SLOWLY INCREASE SX
################################################
mod_map -c run.inpt under_relaxation
mod_map run.inpt dtinv "0.0"
mod_map run.inpt ntstep "1"
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

cd ..

