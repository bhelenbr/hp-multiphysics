#!/bin/bash
#  Runs tests of multigrid convergence for various flow conditions
#

EXECNAME="tri_hp"

# Turn valgrind on/off
let USE_VALGRIND=0  

# Uncomment to set valgrind debugging parameters

VALGRIND_FLAGS+=" --track-origins=yes"
VALGRIND_FLAGS+=" --leak-check=full"
VALGRIND_FLAGS+=" --dsymutil=yes"
VALGRIND="valgrind  ${VALGRIND_FLAGS}"

#Set executable command
if [ "${USE_VALGRIND}" -eq 0 ]; then
	HP="${EXECNAME} run.inpt"
else
	HP="${EXECENV} ${VALGRIND} ${MYPATH}/${EXECNAME}"
fi

cd "$(dirname "$0")"
# Define location of executables
BINDIR=${PWD%/*/*/*/*}/bin
echo ${BINDIR}
export PATH=${PATH}:${BINDIR}


#CASE="FD"
CASE="MG"
#CASE="NO_MG"
#CASE="PMG"

# Make Results directory
if [ -e Results_${CASE} ]; then
	cd Results_${CASE}
else
	mkdir Results_${CASE}
	cd Results_${CASE}
fi
rm -rf *

# copy input files into results directory
cp ../Inputs/* .
cp run_master.inpt run.inpt


case "$CASE" in
FD)
	P="1 2 4"
	GRID="1 2 4 8"
	CFL="0.5 1.0 1.5"
	AR="SQUARE/INOUT_0.1/square SQUARE/INOUT/square SQUARE/INOUT_10/square"
	MGFLAG=1
	INTV="1 20 40 60"
;;

MG)
	mod_map run.inpt ncycle 10
	mod_map run.inpt ntstep 4
	P="1 2 4"
	GRID="2"
	CFL="1.0"
	AR="SQUARE/INOUT/square"
	MGFLAG=1
	INTV="1 2 4 8"
;;
         
NO_MG)
	mod_map run.inpt ncycle 10
	mod_map run.inpt ntstep 4
	P="1 2 4"
	GRID="2"
	CFL="1.0"
	AR="SQUARE/INOUT/square"
	MGFLAG=0
	INTV="1 2 4 8"
;;
PMG)
	mod_map run.inpt ncycle 10
	mod_map run.inpt ntstep 4
	P="2 4"
	GRID="1"
	CFL="1.0"
	AR="SQUARE/INOUT/square"
	MGFLAG=2
	INTV="1 2 4 8"
;;
*)
	echo 'not a valid case';
;;
esac


SUBSP=($P)
SUBSGRID=($GRID)
SUBSCFL=($CFL)
SUBSAR=($AR)
INTERVALS=($INTV)

let node=0
let nar=0
while [ $nar -lt ${#SUBSAR[@]} ]; do
	
	let np=0
	while [ $np -lt ${#SUBSP[@]} ]; do
	   let log2p=${SUBSP[np]}/2
	   mod_map run.inpt log2p ${log2p}

	   let ngrid=0
	   while [ $ngrid -lt ${#SUBSGRID[@]} ]; do
	     let mesh=${SUBSGRID[ngrid]}*4/${SUBSP[np]}
	     mod_map run.inpt b0_mesh ${HOME}/Codes/grids/${SUBSAR[nar]}${mesh}		 

		 if [ $MGFLAG -eq 1 ]; then
			 let lvls=$ngrid-${log2p}+3
			 mod_map run.inpt ngrid ${lvls}
			 mod_map run.inpt extra_finest_levels ${log2p}
		 elif [ $MGFLAG -eq 2 ]; then
		 	 mod_map run.inpt extra_finest_levels ${log2p}
         fi
	   
		 let ncfl=0
		 while [ $ncfl -lt ${#SUBSCFL[@]} ]; do
		 	 CFLS=${SUBSCFL[$ncfl]}
			 mod_map run.inpt C "${CFLS}"

			 
			 name=P${SUBSP[np]}_GRD${SUBSGRID[ngrid]}_AR_${nar}_CFL${ncfl}

			 ##################
			 ### EULER CASE ###
			 ##################
			 if [ -d E_${name} ]; then
			   cd E_${name}
			   rm *
			 else
			   mkdir E_${name};
			   cd E_${name}
			 fi
			 
			 cp ../run.inpt .


			 CMD="${HP} run > ../duh_E_${name}.log"
#			 CMD=qsub ../serial.bash
			 eval $CMD 
			 cd ..;

			 CNVG1=$(grep -v \# duh_E_${name}.log | head -${INTERVALS[0]} | cut -d \  -f 3)
			 CNVG2=$(grep -v \# duh_E_${name}.log | head -${INTERVALS[1]} | tail -1 | cut -d \  -f 3)
			 CNVG3=$(grep -v \# duh_E_${name}.log | head -${INTERVALS[2]} | tail -1 | cut -d \  -f 3)
			 CNVG4=$(grep -v \# duh_E_${name}.log | head -${INTERVALS[3]} | tail -1 | cut -d \  -f 3)
			 echo "0   0    ${SUBSAR[nar]} ${mesh} ${log2p} ${CFLS} ${INTERVALS[0]} ${CNVG1} ${INTERVALS[1]} ${CNVG2} ${INTERVALS[2]} ${CNVG3} ${INTERVALS[3]} ${CNVG4}" >> data.log
			 rm -rf E_${name}
			
			 #################### 
			 ### VISCOUS CASE ###
			 ####################
			 if [ -d V_${name} ]; then
			   cd V_${name}
			   rm *
			 else
			   mkdir V_${name};
			   cd V_${name}
			 fi

			 cp ../run.inpt .
			 mod_map run.inpt b0_mu 100.0
			 CMD="${HP} run > ../duh_V_${name}.log"
#				 CMD=qsub ../serial.bash
			 eval $CMD 
			 cd ..;
			 CNVG1=$(grep -v \# duh_V_${name}.log | head -${INTERVALS[0]} | cut -d \  -f 3)
			 CNVG2=$(grep -v \# duh_V_${name}.log | head -${INTERVALS[1]} | tail -1 | cut -d \  -f 3)
			 CNVG3=$(grep -v \# duh_V_${name}.log | head -${INTERVALS[2]} | tail -1 | cut -d \  -f 3)
			 CNVG4=$(grep -v \# duh_V_${name}.log | head -${INTERVALS[3]} | tail -1 | cut -d \  -f 3)
			 echo "100 0    ${SUBSAR[nar]} ${mesh} ${log2p} ${CFLS} ${INTERVALS[0]} ${CNVG1} ${INTERVALS[1]} ${CNVG2} ${INTERVALS[2]} ${CNVG3} ${INTERVALS[3]} ${CNVG4}" >> data.log
	
			 rm -rf V_${name}
			 
			 #################### 
			 ### UNSTEADY CASE ##
			 ####################
			 if [ -d U_${name} ]; then
			   cd U_${name}
			   rm *
			 else
			   mkdir U_${name};
			   cd U_${name}
			 fi
			 
			 cp ../run.inpt .
			 mod_map run.inpt dtinv 3200.0
			 CMD="${HP} run > ../duh_U_${name}.log"
#				 CMD=qsub ../serial.bash
			 eval $CMD 
			 cd ..;
			 CNVG1=$(grep -v \# duh_U_${name}.log | head -${INTERVALS[0]} | cut -d \  -f 3)
			 CNVG2=$(grep -v \# duh_U_${name}.log | head -${INTERVALS[1]} | tail -1 | cut -d \  -f 3)
			 CNVG3=$(grep -v \# duh_U_${name}.log | head -${INTERVALS[2]} | tail -1 | cut -d \  -f 3)
			 CNVG4=$(grep -v \# duh_U_${name}.log | head -${INTERVALS[3]} | tail -1 | cut -d \  -f 3)
			 echo "0   3200 ${SUBSAR[nar]} ${mesh} ${log2p} ${CFLS} ${INTERVALS[0]} ${CNVG1} ${INTERVALS[1]} ${CNVG2} ${INTERVALS[2]} ${CNVG3} ${INTERVALS[3]} ${CNVG4}" >> data.log
	
			 rm -rf U_${name}
	
			let ncfl=$ncfl+1
			done
	     let ngrid=$ngrid+1
	   done
	   let np=$np+1
	done
	let nar=$nar+1
done
cd ..

opendiff Baseline_${CASE} Results_${CASE}