#!/usr/bin/env bash
set -euo pipefail
#set -x
export FI_PROVIDER=tcp

HP="mpiexec -np 1 tri_hp_petsc"
# Testing accuracy for a case with a singular point

cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/Testing/*}/bin
export PATH=${BINDIR}:${PATH}

mkdir -p Results
cd Results
rm -rf ./*

thetas=(
  "5*_pi/25"
  "5*_pi/15"
  "5*_pi/9"
  "5*_pi/8"
  "5*_pi/7"
  "5*_pi/6"
  "5*_pi/5"
  "5*_pi/4"
  "5*_pi/3"
)

for theta in "${thetas[@]}"; do
	safe_theta=${theta//\//_}   # replace / with _
	safe_theta=${safe_theta//\*/}
	mkdir "$safe_theta"
	cd "$safe_theta"
	cp ../../Inputs/* .
    echo "Running theta = $theta"
    ../../Basictest.command "$theta"
    cd ..
done

