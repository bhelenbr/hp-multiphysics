#!/bin/bash

while getopts ":pt" opt; do
	case $opt in
		p)
			 echo "physical partitions"
			 PHYSICAL="TRUE"
       ;;
		t)
       echo "threaded"
			 THREADED="TRUE"
       ;;
		\?)
       echo 'usage: partition.bash -[p,t] <filename.inpt>'
			 exit 1
			 ;;
	esac
done
shift $((OPTIND-1))

echo $1

MESH=$(grep ^[^#]*esh $1 | cut -d: -f 2)
FTYP=$(grep ^[^#]*iletype $1 | cut -d: -f 2)
BLKS=$(grep ^[^#]*block $1 | cut -d: -f 2)
echo "Partitioning ${MESH} of type ${FTYP} into ${BLKS} parts"

cat $1 | grep -v nblock | grep -v filetype | grep -v b0_mesh > partition.inpt

#echo "text \c"
# or
#echo -n "text "

echo "filetype: 4" >> partition.inpt

if [ -n "${THREADED}" ]; then
# Keep number of blocks the same
	echo "nblock: ${BLKS}" >> partition.inpt
else
# Put 1 block per processor
	echo -n "nblock: " >> partition.inpt
	let p=0
	while [ $p -lt $BLKS ]; do
		echo -n "1 " >> partition.inpt
		let p=$p+1
	done
	echo "" >> partition.inpt
fi

if [ -z "${PHYSICAL}" ]; then
# Copy all boundary information to other blocks
	let p=1
	while [ $p -lt $BLKS ]; do
		grep b0 $1 | grep -v filetype | grep -v mesh | sed "s/b0/b$p/g" >> partition.inpt
		let p=$p+1
	done

# partition mesh
	if [ -n "$FTYP" ];
	then
					echo "tet_mesh -i $FTYP -p $MESH $BLKS | grep -v '#creating' >> partition.inpt"
		tet_mesh -i $FTYP -p $MESH $BLKS | grep -v '#creating'  | grep -v "Info:" >> partition.inpt
	else
					echo "tet_mesh -p $MESH $BLKS | grep -v '#creating' >> partition.inpt"
		tet_mesh -p $MESH $BLKS | grep -v '#creating'  | grep -v "Info:" >> partition.inpt
	fi
else
# divide up physical partitions
	echo "tet_mesh -gp $MESH $BLKS | grep -v '#creating' >> partition.inpt"
	tet_mesh -gp $MESH $BLKS | grep -v '#creating' | grep -v "Info:" >> partition.inpt
fi

# Remove unnecessary communication endpoints (ones only on 2 blocks)
# Crazy unix stuff 
# find all _v lines
# find those that are communications
# get the identifier string
# sort them
# count how many duplicates each has
# sort that again???
# find those with two duplicates
# get just the identifier string again
RMV=$(grep '_v' partition.inpt | grep ': comm' | cut -d _ -f2 | sort | uniq -c | sort | grep '  2 v' | cut -d \  -f 5)

cp partition.inpt temp.inpt
for line in ${RMV}; do
#	echo ${line}
	grep -v ${line} temp.inpt > temp1.inpt
	mv temp1.inpt temp.inpt
done
mv temp.inpt partition.inpt

RMV=$(grep '_e' partition.inpt | grep ': comm' | cut -d _ -f2 | sort | uniq -c | sort | grep '  2 e' | cut -d \  -f 5)

cp partition.inpt temp.inpt
for line in ${RMV}; do
#	echo ${line}
	grep -v ${line} temp.inpt > temp1.inpt
	mv temp1.inpt temp.inpt
done
mv temp.inpt partition.inpt

#RMV=$(grep -w 'Info:' partition.inpt )
#cp partition.inpt temp.inpt
#for line in ${RMV}; do
#	echo ${line}
#	grep -v ${line} temp.inpt > temp1.inpt
#	mv temp1.inpt temp.inpt
#done
#mv temp.inpt partition.inpt

#grep -v 'Info' partition.inpt

sed 's/partition$/comm/g' partition.inpt > comm.inpt
