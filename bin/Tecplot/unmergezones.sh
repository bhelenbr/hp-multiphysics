#!/bin/bash
##
# script by Brian Helenbrook | helenbrk@clarkson.edu
##

for file; do
	LOCS=($(grep -n ZONE $file))
	let count=0
	let bgn=1
	let wc=5
	while [ $wc -lt ${#LOCS[@]} ]; do
		let end=${LOCS[wc]%:ZONE}
		let endm=end-1
		let total=end-bgn
		head -n $endm $file | tail -n $total > ${file%.dat}.$count.dat
		let bgn=end
		let count=count+1
		let wc=wc+5
	done
    tail -n +$bgn $file > ${file%.dat}.$count.dat
done