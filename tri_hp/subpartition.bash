#!/bin/bash

# Script to partition a tri_hp multi-block solution

#echo "text \c"
# or
#echo -n "text "

while getopts ":t" opt; do
	case $opt in
		t)
       echo "threaded"
			 THREADED="TRUE"
       ;;
		\?)
       echo 'usage: partition.bash -[t] <filename.inpt>'
			 exit 1
			 ;;
	esac
done
shift $((OPTIND-1))

cat $1 | grep -v nblock: | grep -v filetype: > partition.inpt
sed -E "s/report: 1/report: 0/g" $1 > temp.inpt
mod_map temp.inpt partition $2
BLKS=$(grep -E "^b[0-9]*_mesh:" $1 | wc -l | tr -d ' ')

#mpiexec -np ${BLKS} ~/bin/tri_hp_petsc temp.inpt -stop_for_debugger
mpiexec -np ${BLKS} ~/bin/tri_hp_petsc temp.inpt
rm temp.inpt

for file in partition_b[0-9].bin; do
	tri_mesh -x $file ${file%.bin}.grd
done

let b=0
rm new.inpt
while [ $b -lt $BLKS ]; do
	grep -E "^[^#]" out_b${b}.log | grep -v matches >> new.inpt
	grep -E "^[^#]" out_b${b}.log | grep matches >> matches.inpt
	let b=$b+1
done
rm out_b*.log

grep -E "b[0-9]*_matches" matches.inpt > block_matches.inpt
grep -E "b[0-9]*_[s,v][0-9]*_matches" matches.inpt > boundary_matches.inpt
rm matches.inpt

let TOTAL=$(ls partition_b*.bin | wc -w | tr -d ' ')

while read tgt src
do
  tgt=${tgt%_matches:}
  tgt=$(echo $tgt | sed 's/^b/B/g')
  grep "$src" partition.inpt | sed "s/$src/$tgt/g" > temp.inpt
  cat temp.inpt >> partition.inpt
done < boundary_matches.inpt
rm boundary_matches.inpt

while read tgt src
do
  tgt=${tgt%_matches:}
  tgt=$(echo $tgt | sed 's/^b/K/g')
  grep "$src" partition.inpt | sed "s/$src/$tgt/g" > temp.inpt
  cat temp.inpt >> partition.inpt
done < block_matches.inpt
rm block_matches.inpt

grep -v "K[0-9][0-9]*_[s,v][0-9][0-9]*" partition.inpt | grep -v "b[0-9][0-9]*" | sed 's/K\([0-9][0-9]*\)/b\1/g' | sed 's/B\([0-9][0-9]*\)/b\1/g' > temp.inpt
rm partition.inpt

if [ -n "${THREADED}" ]; then
# Keep number of blocks the same
	echo "nblock: ${TOTAL}" >> temp.inpt
else
# Put 1 block per processor
	echo -n "nblock: " >> temp.inpt
	let p=0
	while [ $p -lt $TOTAL ]; do
		echo -n "1 " >> temp.inpt
		let p=$p+1
	done
	echo "" >> temp.inpt
fi

RSTRT=$(mod_map -e $1 restart)

if [ -z "${RSTRT}" ]; then
	let RSTRT=1;
fi

rename rstrt${RSTRT} backup${RSTRT} rstrt${RSTRT}*
rename partition rstrt${RSTRT} partition*

let b=0
while [ $b -lt $TOTAL ]; do
	mod_map temp.inpt b${b}_mesh rstrt${RSTRT}_b${b}.bin
	let b=$b+1
done
mod_map -d temp.inpt partition
mv temp.inpt partition.inpt
rm temp.inpt.bak
cat new.inpt >> partition.inpt
rm new.inpt




