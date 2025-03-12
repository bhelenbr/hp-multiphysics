#!/bin/bash

# Number of processors
let NP=1

# cd to the directory where the script resides.
cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/*/*/*/*}/bin
export PATH=${PATH}:${BINDIR}

# Uncomment to add petsc debug flags
#PETSC_FLAGS+=" -info -log_summary -ac-log_summary -memory_info -malloc_log -malloc_info -malloc_debug"
#PETSC_FLAGS+=" -fp_trap"
#PETSC_FLAGS+=" -on_error_attach_debugger gdb"
#PETSC_FLAGS+=" -start_in_debugger gdb"
#PETSC_FLAGS+=" -stop_for_debugger"

# Set executable command
HP="mpiexec -np ${NP} tri_hp_petsc" 
HPAXI="mpiexec -np ${NP} tri_hp_axi_petsc" 

# choose a version:

# version=2D1988
# version=2D1988KL
version=2D2006
# version=2D2006NSL

# Make Results directory
if [ -e Results${version} ]; then
	cd Results${version}
else
	mkdir Results${version}
	cd Results${version}
fi
rm *

# copy input files into results directory
cp ../Inputs/* .

# generate mesh 
tri_mesh generate

${HP} run${version}.inpt ${PETSC_FLAGS}

cd ..

# choose a version:

# version=Axi1988
# version=Axi1988KL
version=Axi2006
# version=Axi2006NSL

# Make Results directory
if [ -e Results${version} ]; then
	cd Results${version}
else
	mkdir Results${version}
	cd Results${version}
fi
rm *

# copy input files into results directory
cp ../Inputs/* .

# generate mesh 
tri_mesh generate

${HPAXI} run${version}.inpt ${PETSC_FLAGS}



