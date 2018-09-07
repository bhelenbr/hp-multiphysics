#!/bin/bash
# Run a free-stream flow with an increasing inlet velocity
# Calculates error in pressure gradient as a function of time step
# Use plot.tank to verify that you get 3'rd order accuracy when using ESDIRK4
# Haven't tried Backwards difference schemes in a long time

HP="tri_hp"

cd "$(dirname "$0")"

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi

rm *
cp ../Inputs/* .

#let backdiff=0

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

opendiff Results/ Baseline/
