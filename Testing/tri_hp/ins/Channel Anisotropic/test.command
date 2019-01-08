#!/bin/bash

# This runs a vertical fully developed channel (periodic flow driven by gravity)

cd "$(dirname "$0")"
BINDIR=${PWD%/*/*/*/*}/bin
export PATH=${PATH}:${BINDIR}
HP="tri_hp_petsc"

#PETSC="-info -on_error_attach_debugger -malloc_log -malloc_info -memory_info"
#PETSC="-info -log_summary -fp_trap -on_error_attach_debugger gdb -malloc_log -malloc_info -malloc_debug -memory_info"
#PETSC="-info -log_summary -fp_trap -start_in_debugger gdb -malloc_log -malloc_info -malloc_debug -memory_info"
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

mpiexec -np 1 ${HP} run.inpt ${PETCSC}


cd ..

opendiff Results Baseline
