#!/bin/bash
# Partitions a tri_hp solution

# cd to the directory where the script resides.
# Necessary for platforms where script can be 
# double clicked from gui
cd "$(dirname "$0")"

# Make Results directory
if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm *

# copy input files into results directory
cp ../Inputs/* .

let NPART=5
~/Codes/tri_hp/subpartition.bash run.inpt ${NPART}

let TOTAL=$(mod_map -e partition.inpt nblock | wc -w | tr -d ' ')
mod_map partition.inpt "growth factor" 3

echo 'b2_v9_hp_type: multi_physics_pnt
b2_v9_b5_v9_matching: 2 0 4 1 5 2
b2_v9_b6_v9_matching: 2 0 4 1 5 2
b4_v9_hp_type: multi_physics_pnt
b4_v9_b5_v9_matching: 2 0 4 1 5 2
b4_v9_b6_v9_matching: 2 0 4 1 5 2
b5_v9_hp_type: multi_physics_pnt
b5_v9_b2_v9_matching: 0 2 1 4 2 5
b5_v9_b4_v9_matching: 0 2 1 4 2 5
b6_v9_hp_type: multi_physics_pnt
b6_v9_b2_v9_matching: 0 2 1 4 2 5
b6_v9_b4_v9_matching: 0 2 1 4 2 5
' >> partition.inpt

sed -i bak 's/: partition/: comm/g' partition.inpt 


mpiexec -np ${TOTAL} tri_hp_petsc partition.inpt 

cd ..

opendiff Baseline/ Results/
