#!/usr/bin/env bash
set -euo pipefail
#set -x
export FI_PROVIDER=tcp

HP="mpiexec -np 1 tri_hp_petsc"

cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/Testing/*}/bin
export PATH=${BINDIR}:${PATH}

mkdir -p Results
cd Results
rm -rf ./*


for n in {8,16}; do
    epsil="exp(-${n})/(1-exp(-${n}))"
    safe_eps="bot${n}"

    mkdir -p "$safe_eps"
    cd "$safe_eps" || exit 1

    cp ../../Inputs/* .
    echo "Running eps = $epsil (dir: $safe_eps)"

    set +e
    ../../Basictest.command "$epsil"
    status=$?
    set -e

    cd ..

    if [[ $status -ne 0 ]]; then
        echo "❌ FAILED for eps = $epsil → removing $safe_eps"
        #rm -rf "$safe_eps"
    else
        echo "✅ SUCCESS for eps = $epsil"
    fi
done

cd ..
./make_plot.command