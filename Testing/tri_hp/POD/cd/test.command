#!/bin/bash
# Runs a convection POD case.
# Single block POD test case.

cd "$(dirname "$0")"

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi

DNS=true
POD_GEN=true
POD_SIM=true

cp ../Inputs/* .

HP="tri_hp"

if [ -n "$DNS" ]; then
	if [ -e DNS ]; then
		rm -rf DNS
	fi
	mkdir DNS
	cd DNS

	${HP} ../dns.inpt
	cd ../
fi

if [ -n "$POD_GEN" ]; then
	if [ -e GEN ]; then
		rm -rf GEN
	fi
	mkdir GEN
	cd GEN

	cp ../DNS/rstrt* .
	${HP} ../pod_gen.inpt
	cd ..
fi

if [ -n "$POD_SIM" ]; then
	if [ -e SIM ]; then
		rm -rf SIM
	fi
	mkdir SIM
	cd SIM

	cp ../GEN/mode* .
	cp ../GEN/rstrt1_* .
	cp ../GEN/coeff1_* .
	${HP} ../pod_sim.inpt
	cd ..
fi

cd ..
opendiff Baseline/ Results/