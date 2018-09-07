#!/bin/bash
# Runs a case with a facet

# uncomment full_test to do the whole thing
# otherwise does one time step only
# remember that "OLD_WAY" has to be undefined in nstage
# Also there are a bunch of switches for turning rotation & swapping off/on
# It hasn't worked with ONE_SIDED not set.  not sure about the others

FULL_TEST=1

if [ -z "$NP" ]; then
	let NP=2;
fi

# Turn valgrind on/off
let USE_VALGRIND=0

# Set executable to run
# Choose this for multigrid
#EXECNAME="tri_hp"
# Choose this for petsc
EXECNAME="tri_hp_petsc"

# Uncomment to add petsc debug flags
#PETSC_FLAGS+=" -info -log_summary -ac-log_summary -memory_info -malloc_log -malloc_info -malloc_debug"
#PETSC_FLAGS+=" -fp_trap"
#PETSC_FLAGS+=" -on_error_attach_debugger gdb"
#PETSC_FLAGS+=" -start_in_debugger gdb"
#PETSC_FLAGS+=" -stop_for_debugger"

# Uncomment to set valgrind debugging parameters
VALGRIND_FLAGS+=" --track-origins=yes"
#VALGRIND_FLAGS+=" --leak-check=yes"
VALGRIND_FLAGS+=" --dsymutil=yes"
 
############################################
# This part should "just work" (in theory)
# Skip to end to modify what is actually done
############################################
# For platform specific valgrind stuff
if [ ${NP} -gt 0 ]; then
	if [ -e "${HOME}/Packages/lib/valgrind/libmpiwrap-amd64-darwin.so" ]; then
		export LD_PRELOAD="${HOME}/Packages/lib/valgrind/libmpiwrap-amd64-darwin.so"
		VALGRIND_FLAGS+=" --suppressions=${HOME}/Packages/share/openmpi/openmpi-valgrind.supp"
	else
		VALGRIND_FLAGS+=" --suppressions=/usr/local/src/mvapich2-1.9a/src/mpid/ch3/channels/mrail/src/hwloc/contrib/hwloc-valgrind.supp"
	fi
	EXECENV="mpiexec -np ${NP} ${MF}"
fi

VALGRIND="valgrind  ${VALGRIND_FLAGS}"
MYPATH="${HOME}/bin"

#Set executable command
if [ ${USE_VALGRIND} -eq 0 ]; then
	HP="${EXECENV} ${MYPATH}/${EXECNAME} run.inpt ${PETSC_FLAGS}"
else
	HP="${EXECENV} ${VALGRIND} ${MYPATH}/${EXECNAME} run.inpt ${PETSC_FLAGS}"
fi

#############################################
# End of automated part HP is now the correct
# string needed to run the executable     
# can just copy and paste above part into shell if you want
# to run the executable without using this script
# Paste above into shell to define HP
# then type ${HP} to run
#############################################


# cd to the directory where the script resides.
# Necessary for platforms where script can be 
# double clicked from gui
cd "$(dirname "$0")"

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
delete_data.bash 2 0


################################################
# SX IS THE PULL SPEED
# NOW SLOWLY INCREASE SX
################################################
mod_map -c run.inpt under_relaxation
mod_map run.inpt dtinv "0.0"
mod_map run.inpt ntstep "1"
SX="$(mod_map -e run.inpt sx | cut -d. -f1)0"
SXMAX=1000000
FACTOR="1.05"
let ATTEMPTS=0
let MAXATTEMPTS=3

let CURVE_START=${RESTART}+1

while [ ${ATTEMPTS} -lt ${MAXATTEMPTS} ] && [ $(echo "${SX} < ${SXMAX}" | bc -l) -eq 1 ]; do
	SXOLD=${SX}
	SX=$(echo "scale=6;${SX}*${FACTOR}/1" | bc -l)
	echo "SX ${SX}"
	
	mod_map run.inpt sx "${SX}e-7/(d0*tsi)"
	mod_map run.inpt restart ${RESTART}

	${HP}
	if [ "$?" -ne "0" ]; then
	  rm core*
	  rm neg*
	  rm abort*
	  SX=${SXOLD}
	  FACTOR=$(echo "scale=5;(${FACTOR}-1)/2 +1" | bc -l)
	  let ATTEMPTS=${ATTEMPTS}+1
	  let RESTART=${RESTART}-1
	else
		# TEST TO MAKE SURE CONVERGED
		RES=$(grep ^[1-9] out_b0.log | tail -1 | cut -d\  -f 3)
		RES=${RES/[eE]/*10^\(}\)
		RES=${RES/\([+]/\(}
		CNVG=$(echo "${RES} < 1.0*10^(-5)" | bc -l)
		if [ "$CNVG" != "1" ]; then
	  		SX=${SXOLD}
	  		FACTOR=$(echo "scale=5;(${FACTOR}-1)/2 +1" | bc -l)
	  		let ATTEMPTS=${ATTEMPTS}+1
	  		let RESTART=${RESTART}-1
		fi
	fi
	let RESTART=${RESTART}+1
done
rm core*
rm abort*
rm neg*
rm rstrt*.bin

cd ..

