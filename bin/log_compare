#!/bin/bash -f
let i=1
let ip=$i+1
while [ -e out_b0.${ip}.log ]; do
	file1=out_b0.${i}.log
	file2=out_b0.${ip}.log
	LINES1=$(grep -n -m 1 TIMESTEP ${file1} | cut -f 1 -d:)
	LINES2=$(grep -n -m 1 TIMESTEP ${file2} | cut -f 1 -d:)
	let LINES1=${LINES1}-1
	let LINES2=${LINES2}-1
	head -${LINES1} ${file1} | sort | uniq > tmp1.dat
	head -${LINES2} ${file2} | sort | uniq > tmp2.dat
	echo ${file2}  >> change.log
	diff -w  tmp1.dat tmp2.dat >> change.log
	tail -1 ${file2} >> change.log
    let i=$i+1
    let ip=$i+1
done
file1=out_b0.${i}.log
file2=out_b0.log
LINES1=$(grep -n -m 1 TIMESTEP ${file1} | cut -f 1 -d:)
LINES2=$(grep -n -m 1 TIMESTEP ${file2} | cut -f 1 -d:)
let LINES1=${LINES1}-1
let LINES2=${LINES2}-1
head -${LINES1} ${file1} | sort | uniq > tmp1.dat
head -${LINES2} ${file2} | sort | uniq > tmp2.dat
echo ${file2}  >> change.log
diff -w  tmp1.dat tmp2.dat >> change.log
tail -1 ${file2} >> change.log
let i=$i+1
let ip=$i+1





# LINES1=$(grep -n -m 1 TIMESTEP $1 | cut -f 1 -d:)
# LINES2=$(grep -n -m 1 TIMESTEP $2 | cut -f 1 -d:)
# 
# let LINES1=${LINES1}-1
# let LINES2=${LINES2}-1
# head -${LINES1} $1 | sort | uniq > tmp1.dat
# head -${LINES2} $2 | sort | uniq > tmp2.dat
# diff tmp1.dat tmp2.dat | sort | uniq > $1_$2.changes
# tail -1 $2 >> $1_$2.changes

#opendiff tmp1.dat tmp2.dat