#!/bin/bash
# This runs the Taylor Green vortex problem (double periodic)
# and calculates the error using different time steps
# Only 5 timesteps are needed to hit the maximum error
# even as the time step is decreased
# output is file cnvg.dat which can be plotted using
# convergence plots.tank
# Can also get contour plots of the error using
# plot.tank verify that you get third order accuracy using ESDIRK4
# pressure doesn't converge because of offset by a constant
# Last check gave convergence fit of y = 0.0452535 x^{2.76978}

HP="tri_hp"

cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/*/*/*/*/*}/bin
export PATH=${PATH}:${BINDIR}

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi

rm *
cp ../Inputs/* .

# let backdiff=1

let dtinv=30
let ntstep=5
let dtc=1

# if [ -d ${backdiff} ]; then
#   let dtinv=$dtinv*3
#   let ntstep=$ntstep*3
# fi
	
while [ $dtc -lt 5 ]; do
	mod_map run.inpt logfile output${dtc}
	mod_map run.inpt dtinv ${dtinv}
	mod_map run.inpt ntstep ${ntstep}
	${HP} run
	
	grep L_2 output${dtc}_b0.log | sed "s/\#L_2://g" | sed "s/L_inf \([-+.e0-9]*\)[ ]*[0-9]*/\1/g" > err${dtc}
    l2umax=$(cut -f 2 -d" " err${dtc} | awk -v max=0 '{if($1>max){max=$1}}END{print max}')
    liumax=$(cut -f 3 -d" " err${dtc} | awk -v max=0 '{if($1>max){max=$1}}END{print max}')
    l2vmax=$(cut -f 4 -d" " err${dtc} | awk -v max=0 '{if($1>max){max=$1}}END{print max}')
    livmax=$(cut -f 5 -d" " err${dtc} | awk -v max=0 '{if($1>max){max=$1}}END{print max}')
    l2pmax=$(cut -f 6 -d" " err${dtc} | awk -v max=0 '{if($1>max){max=$1}}END{print max}')
    lipmax=$(cut -f 7 -d" " err${dtc} | awk -v max=0 '{if($1>max){max=$1}}END{print max}')
    echo "$l2umax $liumax $l2vmax $livmax $l2pmax $lipmax" >> cnvg.dat 
    
    
	let dtc=${dtc}+1
	let dtinv=2*${dtinv}
	#let ntstep=2*${ntstep}
done

cd ..

opendiff Results/ BASELINE/
