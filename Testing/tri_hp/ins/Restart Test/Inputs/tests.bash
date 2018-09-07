SUBSP=($P)

	
let np=0
while [ $np -lt ${#SUBSP[@]} ]; do
	let log2p=${SUBSP[np]}/2
	mod_map run.inpt log2p ${log2p}

	if [ $MGFLAG -eq 1 ]; then
	 mod_map run.inpt extra_finest_levels ${log2p}
	elif [ $MGFLAG -eq 2 ]; then
	 mod_map run.inpt extra_finest_levels ${log2p}
	fi

	name=P${SUBSP[np]}

	#################### 
	### STEADY CASE ###
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
	CMD="${HP} run > duh_V_${name}.log"
	eval $CMD 
	NLINES=$(grep -n -m 1 'TIMESTEP: 4' duh_V_${name}.log | cut -d: -f1)
	tail -n +${NLINES} duh_V_${name}.log | head -3 > duh_V_${name}_last.log

	mod_map run.inpt restart 3
	mod_map run.inpt ntstep 1
	CMD="${HP} run > duh_V_${name}_rstrt.log"
	eval $CMD
	NLINES=$(grep -n -m 1 'TIMESTEP: 4' duh_V_${name}_rstrt.log | cut -d: -f1)
	tail -n +${NLINES} duh_V_${name}_rstrt.log | head -3 > duh_V_${name}_rstrt_last.log

	diff duh_V_${name}_last.log duh_V_${name}_rstrt_last.log > ../duh_V_${name}_diff.log
	rm duh_V_${name}_last.log duh_V_${name}_rstrt_last.log
	cd ..;

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
	CMD="${HP} run > duh_U_${name}.log"
	eval $CMD 
	NLINES=$(grep -n -m 1 'TIMESTEP: 4' duh_U_${name}.log | cut -d: -f1)
	tail -n +${NLINES} duh_U_${name}.log | head -3 > duh_U_${name}_last.log


	mod_map run.inpt restart 3
	mod_map run.inpt ntstep 1
	CMD="${HP} run > duh_U_${name}_rstrt.log"
	eval $CMD
	NLINES=$(grep -n -m 1 'TIMESTEP: 4' duh_U_${name}_rstrt.log | cut -d: -f1)
	tail -n +${NLINES} duh_U_${name}_rstrt.log | head -3 > duh_U_${name}_rstrt_last.log

	diff duh_U_${name}_last.log duh_U_${name}_rstrt_last.log > ../duh_U_${name}_diff.log
	rm duh_U_${name}_last.log duh_U_${name}_rstrt_last.log
	cd ..;	

		
	let np=$np+1
done
cd ..
