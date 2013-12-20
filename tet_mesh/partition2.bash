#!/bin/bash

#  partition2.bash
#  tet_mesh
#
#  Created by Brian Helenbrook on 12/18/13.
#

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

MESH=$(grep ^b0_mesh $1 | cut -d: -f 2)
FTYP=$(grep ^filetype $1 | cut -d: -f 2)
BLKS=$(grep ^nblock $1 | cut -d: -f 2)
echo "Partitioning ${MESH} of type ${FTYP} into ${BLKS} parts"

cat $1 | grep -v nblock | grep -v filetype | grep -v b0_ > partition.inpt

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
# partition mesh
	if [ -n "$FTYP" ];
	then
					echo "tet_mesh -i $FTYP -p $MESH $BLKS | grep -v Info > boundary_mapping.txt"
		tet_mesh -i $FTYP -p $MESH $BLKS | grep -v "Info"  > boundary_mapping.txt
	else
					echo "tet_mesh -p $MESH $BLKS | grep -v Info > boundary_mapping.txt"
		tet_mesh -p $MESH $BLKS | grep -v "Info"  > boundary_mapping.txt
	fi

	# get the lines defining the new partition boundaries
	grep partition boundary_mapping.txt >> partition.inpt

	# parse boundary mapping file and apply changes to create partition.inpt
	OIFS=${IFS}
	IFS=$'\n'
	RMV=$(grep -v Info boundary_mapping.txt | grep -v partition)
	for line in ${RMV}; do
		NEW=$(echo ${line} | cut -d \  -f 1)
		OLD=$(echo ${line} | cut -d \  -f 2)
		grep ${OLD} $1 | sed "s/${OLD}/${NEW}/g" >> partition.inpt
	done
	IFS=${OIFS}
	echo "b0_mesh: partition_b0" >> partition.inpt

	# Copy all boundary information to other blocks
	let p=1
	while [ $p -lt $BLKS ]; do
		grep b0 partition.inpt | sed "s/b0/b$p/g" >> partition.inpt
		let p=$p+1
	done

else
# divide up physical partitions
	echo "tet_mesh -gp $MESH $BLKS | grep -v '#creating' >> partition.inpt"
	tet_mesh -gp $MESH $BLKS | grep -v '#creating' | grep -v "Info:" >> partition.inpt
fi
