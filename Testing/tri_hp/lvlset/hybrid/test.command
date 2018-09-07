#!/bin/bash
# Runs a grid convergence study with a hybrid_pt at the inflow
# and outflow periodic boundaries

# To tun this case, must make sure flow is not updated (u=1)
# Have to brute force make sure in both level-set and ins in tstep.cpp: setup_preconditioner
# Also turn on #define l2error in bdry_ins.h

# Other permutations are to muck with the level-set boundary condition
# can uncomment return in hybrid::update to simply freeze values on incoming boundaries
# can set reinit iteration to in input file to turn on level-set re-initialization

# To run test case comment 
CNVG=true

if [ -n "$CNVG" ]; then
	NGRID=4
else
	NGRID=2
fi


cd "$(dirname "$0")"

if [ -e $HOME/Codes/tri_hp/build/Release/tri_hp ]; then
  HP="$HOME/Codes/tri_hp/build/Release/tri_hp"
else
  HP="$HOME/bin/tri_hp"
fi

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm -rf *

cp ../Inputs/* .

# Make the mesh
cp run.inpt generate.inpt
mod_map generate.inpt ntstep 2
mod_map generate.inpt adapt 1
mod_map generate.inpt b0_filetype 8
mod_map generate.inpt b1_filetype 8
mod_map generate.inpt b2_filetype 8
mod_map generate.inpt b0_growth\ factor 1000.0
mod_map generate.inpt b1_growth\ factor 1000.0
mod_map generate.inpt b2_growth\ factor 1000.0
mod_map generate.inpt b0_s5_type symbolic_comm
mod_map generate.inpt b1_s5_type symbolic_comm
mod_map generate.inpt ngrid 1
mod_map generate.inpt log2p 0
mod_map generate.inpt extra_finest_levels 0
tri_mesh generate
cp data2_b0.grd botleft1.grd
cp data2_b1.grd topleft1.grd
cp data2_b2.grd right1.grd
rm data*.grd

let ngrid=1
let dtinv=4
while [ $ngrid -le ${NGRID} ]; do
	let dtinv=${dtinv}*2
	let ntstep=4*dtinv
	mkdir grid${ngrid}
	cd grid${ngrid}
	cp ../run.inpt .
	cp ../*.grd .
	mod_map run.inpt ngrid ${ngrid}
	mod_map run.inpt dtinv ${dtinv}
	if [ -n "$CNVG" ]; then
		mod_map run.inpt ntstep ${ntstep}
	fi
	${HP} run
	grep Contour out_b2.log | cut -d\    -f3 > l2_zero${ngrid}.dat
	grep b0_s5\ l2 out_b0.log | cut -d\    -f4 > l2_iface${ngrid}.dat
	mergelines l2_zero${ngrid}.dat l2_iface${ngrid}.dat > l2error${ngrid}.dat
	rm l2_zero${ngrid}.dat l2_iface${ngrid}.dat
	mv l2error${ngrid}.dat ..
	let ngp=${ngrid}+1
	cd ..
	tri_mesh -r botleft1 botleft1
	tri_mesh -r topleft1 topleft1
	tri_mesh -r right1 right1
	let ngrid=${ngp}
done

cd ..
