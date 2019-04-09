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



# To extract debugging statements printed in residual with the keyword "komega"
# grep komega out_b0.log | sed 's/komega//' | sort -g | uniq >  komega.log

# This is for when rsdl_debug is set to 1.  It will extract the residuals and merge them with a data0_b0.dat file
# as extra columns so that you can plot what the residuals look like.  I still had to remove extra spaces from the file
# to get it to work
# echo ' ' > temp.dat
# grep 'b0 v:' out_b0.log | sed 's/b0 v://' >> temp.dat
# paste  data0_b0.dat temp.dat | sed "s/^$(printf '\t')$//g" > residual_b0.dat

cd ..

opendiff Results Baseline
