#!/bin/bash

MESH=$(grep ^[^#]*esh $1 | cut -d: -f 2)
FTYP=$(grep ^[^#]*iletype $1 | cut -d: -f 2)
BLKS=$(grep ^[^#]*block $1 | cut -d: -f 2)
echo "Partitioning ${MESH} of type ${FTYP} into ${BLKS} parts"

cat $1 | grep -v nblock | grep -v filetype | grep -v b0_mesh > partition.inpt

#echo "text \c"
# or
#echo -n "text "

echo "filetype: 3" >> partition.inpt
echo -n "nblock: " >> partition.inpt
let p=0
while [ $p -lt $BLKS ]; do
	echo -n "1 " >> partition.inpt
	let p=$p+1
done
echo "" >> partition.inpt

let p=1
while [ $p -lt $BLKS ]; do
	grep b0 $1 | grep -v filetype | grep -v mesh | sed "s/b0/b$p/g" >> partition.inpt
	let p=$p+1
done

#echo $FTYP
#echo $MESH
#echo $BLKS
	
if [ -n "$FTYP" ]; 
then
        echo "tri_mesh -i $FTYP -p $MESH $BLKS | grep -v '#creating' >> partition.inpt"
	tri_mesh -i $FTYP -p $MESH $BLKS | grep -v '#creating' >> partition.inpt
else
        echo "tri_mesh -p $MESH $BLKS | grep -v '#creating' >> partition.inpt"
	tri_mesh -p $MESH $BLKS | grep -v '#creating' >> partition.inpt
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
sed 's/partition$/comm/g' partition.inpt > comm.inpt
