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
# check rates.dat Error.pdf
# Need to recompile to run BD schemes

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
rm -rf *

# let backdiff=1

let dtinv=30
let ntstep=5

# if [ -d ${backdiff} ]; then
#   let dtinv=$dtinv*3
#   let ntstep=$ntstep*3
# fi

names=("DIRK1" "DIRK2" "DIRK3" "DIRK4" "AM1" "BD1" "BD2" "BD3")
SCHEMES=(6 7 8)
SCHEMES=(1 2 3 4 5)


let nscheme=0
while [ $nscheme -lt ${#SCHEMES[@]} ]; do
	let sch=${SCHEMES[nscheme]}
	name=${names[sch-1]}
	mkdir ${name}
	cd ${name}
	cp ../../Inputs/* .
	mod_map run.inpt time_scheme ${sch}
	
	let dtinv=30
	let ntstep=5
	let dtc=1
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
		echo "$l2umax $liumax $l2vmax $livmax" >> cnvg.dat 
		
		
		let dtc=${dtc}+1
		let dtinv=2*${dtinv}
		#let ntstep=2*${ntstep}
	done
	
	cd ..
	
	../make_plot.command ${name} > ${name}/rates.dat
	let nscheme=${nscheme}+1
done

cd ..

opendiff Results/ Baseline/
