# Tests GCL condition for curved interior edges
# In fact, this test fails. GCL is only satisfied
# For linear elements

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

# Using Petsc parallel
HP="mpiexec -np 2 tri_hp_petsc run.inpt"

#PETSC="-info -on_error_attach_debugger -malloc_log -malloc_info -memory_info"
#PETSC="-info -log_summary -fp_trap -on_error_attach_debugger gdb -malloc_log -malloc_info -malloc_debug -memory_info"
#PETSC="-info -log_summary -fp_trap -start_in_debugger gdb -malloc_log -malloc_info -malloc_debug -memory_info"
#PETSC="-stop_for_debugger"

cp ../Inputs/* .

${HP} ${PETSC}

cd ..
opendiff Results Baseline
