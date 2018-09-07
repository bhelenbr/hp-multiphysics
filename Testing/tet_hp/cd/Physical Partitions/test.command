#!/bin/bash

cd "$(dirname "$0")"

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm *

# Solve for the modes
cp ../Inputs/* .

${HOME}/Codes/tet_mesh/partition.bash -p -t run.inpt
#${HOME}/Codes/tet_mesh/partition2.bash run.inpt


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

CMMD="${EXECENV} ${MYPATH}/${EXECNAME} comm.inpt ${FLAGS}"

${CMMD}

cd ../
opendiff Baseline/ Results/
