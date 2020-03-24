#!/bin/bash

# 2D turbulent planar stagnation flow used to test k-ln(omg)1988-KL model with turbulence sustaining terms turned on

# Number of processors
let NP=1

# cd to the directory where the script resides.
cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/*/*/*/*/*}/bin
export PATH=${PATH}:${BINDIR}

# Uncomment to add petsc debug flags
#PETSC_FLAGS+=" -info -log_summary -ac-log_summary -memory_info -malloc_log -malloc_info -malloc_debug"
#PETSC_FLAGS+=" -fp_trap"
#PETSC_FLAGS+=" -on_error_attach_debugger gdb"
#PETSC_FLAGS+=" -start_in_debugger gdb"
#PETSC_FLAGS+=" -stop_for_debugger"

# Set executable command
HP="mpiexec -np ${NP} tri_hp_axi_petsc run.inpt ${PETSC_FLAGS}"

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

# generate mesh 
tri_mesh generate

# Run the executable 
${HP}




