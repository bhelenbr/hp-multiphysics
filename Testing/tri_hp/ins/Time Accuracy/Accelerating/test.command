#!/bin/bash
# Run a free-stream flow with an increasing inlet velocity
# Calculates error in pressure gradient as a function of time step
# Haven't tried Backwards difference schemes in a long time
# Check rates.dat & Error.pdf for convergence

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

cp ../Inputs/*.grd .
# Make meshes
tri_mesh -m 'x0+0.01*x1,x1' square1.grd distorted1.grd

let ngrid=1
while [ $ngrid -lt 3 ]; do
	let ngp=${ngrid}+1
	tri_mesh -r distorted${ngrid}.grd distorted${ngp}.grd
	let ngrid=${ngrid}+1
done

let ngrid=1
while [ $ngrid -lt 4 ]; do
	tri_mesh -m 'x0-0.01*x1,x1' distorted${ngrid}.grd square${ngrid}.grd
	let ngrid=${ngrid}+1
done
rm distorted*


let backdiff=0

let dtinv=20
let ntstep=4
let dtc=1

if [ $backdiff -gt 0 ]; then
  let dtinv=$dtinv*3
  let ntstep=$ntstep*3
fi
	
while [ $dtc -lt 5 ]; do
	mod_map run.inpt logfile output${dtc}
	mod_map run.inpt dtinv ${dtinv}
	mod_map run.inpt ntstep ${ntstep}
	${HP} run

	grep L_2 output${dtc}_b0.log | tail -1 | cut -d\  -f2,4,6,8,10,12 >> cnvg.dat

	let dtc=${dtc}+1
	let dtinv=2*${dtinv}
	# let ntstep=2*${ntstep}
done

cd ..

./make_plot.command > Results/rates.dat

opendiff Results/ Baseline/
