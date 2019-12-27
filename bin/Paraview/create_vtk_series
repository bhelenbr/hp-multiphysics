#!/bin/zsh 

usage() {                                      # Function: Print a help message.
  echo 'usage: create_series.vtk -b <nblock> -n <start> <end> <interval>'  1>&2 
  echo 'or:    create_series.vtk -b <nblock> <b0 filelist>' 1>&2 
  echo '-o <name> will cause the output files to be named name#.vtk.series'
  echo '-p <name> will look for files of type name#_b#.vtk instead of data_#_b#.vtk'
}

exit_abnormal() {                              # Function: Exit with error.
  usage
  exit 1
}

let NBLOCKS=1
FNAME='data'
ONAME='b'

while getopts "hnb:o:p:" opt; do
   case $opt in
   n)
   	let SEQUENCE=1
   	;;
   b)
   	NBLOCKS=${OPTARG}                        # Set $NBLOCKS to specified value.
    re_isanum='^[0-9]+$'                     # Regex: match whole numbers only
	if ! [[ $NBLOCKS =~ $re_isanum ]] ; then  # if $NBLOCKS not a whole number:
		echo "Error: NBLOCKS must be a positive, whole number."
		exit_abnormally
		exit 1
	elif [ $NBLOCKS -eq "0" ]; then            # If it's zero:
		echo "Error: NBLOCKS must be greater than zero."
		exit_abnormal                          # Exit abnormally.
	fi
	;;
	o)
	ONAME=${OPTARG}
	;;
	p)
	FNAME=${OPTARG}
	;;
   \?) 
          exit_abnormal
   esac
done
shift $((OPTIND-1))

echo ${ONAME}
echo ${FNAME}


if [ "${SEQUENCE}" -eq 1 ]; then
	# Check interval
	if [ $# -eq 3 ]; then
		let INT=$3
	else
		let INT=1
	fi

	for (( b = 0; b<$NBLOCKS; ++b )) do
		echo '{
		  "file-series-version" : "1.0",
		  "files" : [' > ${ONAME}${b}.vtk.series
	
		for (( c = $1; c< $2; c+=${INT} )) do
			echo "{ \"name\" : \"${FNAME}${c}_b${b}.vtk\", \"time\" : $c },"  >> ${ONAME}${b}.vtk.series
		done
		echo "{ \"name\" : \"${FNAME}${2}_b${b}.vtk\", \"time\" : $c }"  >> ${ONAME}${b}.vtk.series
		echo "]
		}"  >> ${ONAME}${b}.vtk.series
	done
else
	for (( b = 0; b<$NBLOCKS; ++b )) do
		echo '{
		  "file-series-version" : "1.0",
		  "files" : [' > ${ONAME}${b}.vtk.series
	
		let NFILES=$#
		SORTED_FILES=$(for file in "$@"; do
    		echo $file
		done | sort -V | tr '\n' ' ')
		IFS=' ' read -r -A array <<< "$SORTED_FILES"
		for (( c = 1; c< ${NFILES}; c+=1 )) do
			file=${array[$c]%_b*.vtk}
			echo "{ \"name\" : \"${file}_b${b}.vtk\", \"time\" : ${file#${FNAME}} },"  >> ${ONAME}${b}.vtk.series
		done
		file=${array[$NFILES]%_b*.vtk}
		echo "{ \"name\" : \"${file}_b${b}.vtk\", \"time\" : ${file#${FNAME}} }"  >> ${ONAME}${b}.vtk.series
		echo "]
		}"  >> ${ONAME}${b}.vtk.series
	done
fi