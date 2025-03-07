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
# Define location of executables
BINDIR=${PWD%/*/*/*/*}/bin
echo ${BINDIR}
export PATH=${PATH}:${BINDIR}

HP="mpiexec -np 1 tri_hp_petsc"
export HP

P="1 2 4"
GRID="4"
CFL="1.0"
AR="square"
MGFLAG=0
INTV="2 3 4 5"

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

cp ../Inputs/* .
# Make meshes
tri_mesh -m 'x0+0.01*x1,x1' square1.grd distorted1.grd

let ngrid=1
while [ $ngrid -lt 4 ]; do
	let ngp=${ngrid}+1
	tri_mesh -r distorted${ngrid}.grd distorted${ngp}.grd
	let ngrid=${ngrid}+1
done

let ngrid=1
while [ $ngrid -lt 5 ]; do
	tri_mesh -m 'x0-0.01*x1,x1' distorted${ngrid}.grd square${ngrid}.grd
	tri_mesh -m 'x0,x1*0.1' square${ngrid}.grd narrow${ngrid}.grd
	tri_mesh -m 'x0,x1*10.0' square${ngrid}.grd wide${ngrid}.grd
	let ngrid=${ngrid}+1
done


# Test just flow convergence
mkdir Flow
cd Flow
cp ../* .
./tests.bash
cd ../

# Test deforming mesh convergence
mkdir Deforming
cd Deforming
cp ../* .
mod_map -u run.inpt mesh_movement
./tests.bash
cd ../

# Test periodic mesh
mkdir Periodic
cd Periodic
cp ../* .
mod_map run.inpt b0_s1_type prdc
mod_map -c run.inpt b0_s1_hp_type
./tests.bash
cd ../

cd ..
opendiff Results Baseline