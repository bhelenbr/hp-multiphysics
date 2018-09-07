#!/bin/bash
# Runs a cubic solution with a source
# Should get analytic solution I think
# Need to double check this

cd "$(dirname "$0")"

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm *
cp ../Inputs/* .

# For petsc
#EXECNAME="tet_hp_petsc"
#EXECENV="${HOME}/Packages/bin/mpiexec -np 2"
#FLAGS="-info -log_summary -fp_trap -on_error_attach_debugger gdb -malloc_log -malloc_info -malloc_debug -memory_info"
#FLAGS="-start_in_debugger gdb"
#FLAGS="-stop_for_debugger"

# For multigrid
EXECNAME="tet_hp"

OSXPATH="${HOME}/Codes/tet_hp/DerivedData/tet_hp/Build/Products/"
if [ -e $OSXPATH ]; then
	if [ ${OSXPATH}/Debug/${EXECNAME} -nt ${OSXPATH}/Release/${EXECNAME} ]; then
		MYPATH=${OSXPATH}/Debug/
	else
		MYPATH=${OSXPATH}/Release/
	fi
else
  MYPATH="$HOME/bin/"
fi

CMMD="${EXECENV} ${MYPATH}/${EXECNAME} run ${FLAGS}"

${CMMD}

cd ../
opendiff Baseline/ Results/
