#!/bin/bash 

for file in $@; do
	# Total lines in file
	TLINES=$(wc -l ${file} | cut -d\  -f 3)
	echo ${TLINES}
	
	# Line at which 1D zone starts
	LINE=$(grep -n ZONE $file | head -2 | tail -1 | cut -d: -f1)
	let LINE=${LINE}-1
	echo ${LINE}
	head -${LINE} ${file} > ${file%.dat}_cut.dat

	let NTAIL=${TLINES}-${LINE}
	echo ${NTAIL}
	tail -${NTAIL} ${file} > ${file%.dat}_1d.dat
done

