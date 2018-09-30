# Tests periodic boundary conditions with a 4 way vertex point
# Answer should be a vertical free-stream

# This converges, but there seems to be some error in the Jacobian
# because it doesn't seem to converge quadratically
# (at least it works though, sort of a tough case)

cd "$(dirname "$0")"

# Define location of executables
BINDIR=${PWD%/*/*/*/*}/bin
echo ${BINDIR}
export PATH=${PATH}:${BINDIR}

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm -rf *

cp ../Inputs/* .

tri_mesh generate.inpt
rm data*.grd

# CASE 2 DEBUG TEST WITH PETSC MULTIGRID
#LOG2P="0"
LOG2P="0 1 2"
LOG2P=(${LOG2P})

# Using Petsc parallel
HP="mpiexec -np 2 tri_hp_petsc run.inpt"

#PETSC="-info -on_error_attach_debugger -malloc_log -malloc_info -memory_info"
#PETSC="-info -log_summary -fp_trap -on_error_attach_debugger gdb -malloc_log -malloc_info -malloc_debug -memory_info"
#PETSC="-info -log_summary -fp_trap -start_in_debugger gdb -malloc_log -malloc_info -malloc_debug -memory_info"
#PETSC="-stop_for_debugger"

# Test just flow convergence
mkdir Flow
cd Flow
cp ../../Inputs/* .
let np=0
while [ $np -lt ${#LOG2P[@]} ]; do
   log2p=${LOG2P[np]}
   mod_map run.inpt log2p ${log2p}
   name="log2p=${log2p}"
   mod_map run.inpt logfile ${name}
   ${HP} ${PETSC}
   let np=${np}+1
done
cd ../

# Test deforming mesh convergence
mkdir Deforming
cd Deforming
cp ../../Inputs/* .
mod_map -u run.inpt mesh_movement
let np=0
while [ $np -lt ${#LOG2P[@]} ]; do
   log2p=${LOG2P[np]}
   mod_map run.inpt log2p ${log2p}
   name="log2p=${log2p}"
   mod_map run.inpt logfile ${name}
   ${HP} ${PETSC}
   let np=${np}+1
done
cd ../

# Test periodic mesh
mkdir Periodic
cd Periodic
cp ../../Inputs/* .
mod_map run.inpt b0_s2_type prdc
mod_map -c run.inpt b0_s2_hp_type
mod_map run.inpt b1_s4_type prdc
mod_map -c run.inpt b1_s4_hp_type
mod_map -u run.inpt b0_v1_type
mod_map -u run.inpt b1_v1_type
mod_map -u run.inpt b0_s2_type
mod_map -u run.inpt b1_s4_type
let np=0
while [ $np -lt ${#LOG2P[@]} ]; do
   log2p=${LOG2P[np]}
   mod_map run.inpt log2p ${log2p}
   name="log2p=${log2p}"
   mod_map run.inpt logfile ${name}
   ${HP} ${PETSC}
   let np=${np}+1
done
cd ../

cd ..
opendiff Results Baseline
