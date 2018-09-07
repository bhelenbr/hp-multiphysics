# This tests petsc convergence
# for different flow conditions
# and for 
# 1) flow 
# 2) deforming mesh & flow 
# 3) periodic single block b.c.'s
# For periodic meshes can't use resolution less than 4
# because elements touch across the domain messing up the
# calculation of the number of non-sparse points.
# Baseline_OLD got slightly better convergence (GRID="2")
# Probably because of changes to numerical evaluation of Jacobian

cd "$(dirname "$0")"

HP="mpiexec -np 1 tri_hp_petsc"
export HP

P="1 2 4"
GRID="4"
CFL="1.0"
AR="SQUARE/INOUT/square"
MGFLAG=0
INTV="1 2 4 8"

export P
export GRID
export CFL
export AR
export MGFLAG
export INTV

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm -rf *

# Test just flow convergence
mkdir Flow
cd Flow
cp ../../Inputs/* .
./tests.bash
cd ../

# Test deforming mesh convergence
mkdir Deforming
cd Deforming
cp ../../Inputs/* .
mod_map -u run.inpt mesh_movement
./tests.bash
cd ../

# Test periodic mesh
mkdir Periodic
cd Periodic
cp ../../Inputs/* .
mod_map run.inpt b0_s1_type prdc
mod_map -c run.inpt b0_s1_hp_type
./tests.bash
cd ../

cd ..
opendiff Results Baseline