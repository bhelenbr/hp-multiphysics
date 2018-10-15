#!/bin/bash

while getopts ":drst" opt; do
   case $opt in
      d) Command="Delete";;
      r) Command="Revert";;
      s) Command="Search";;
      t) Command="Test";;

      \?) echo 'usage: global_inpt_change.command -[drst] <Directory> <Argument1> <Argument2>'
      	  echo 'd: delete lines matching argument 1'
      	  echo 'r: Revert change argument is time stamp'
      	  echo 's: Search for argument 1 and list file names'
      	  echo 't: Test change, arguments 1/2 are search/replace strings'
      	  echo 'default is to do find/replace'
          exit 1
   esac
done
shift $((OPTIND-1))

timestamp=`date +%Y%m%d%H%M`

# Using ;'s as delimiter for replace string
# because it is already a special character anyway
# avoids problems with replace strings with a /
case "${Command}" in 
	Delete)
		find $1 -name '*.inpt' | while read file
		do
  			if grep "$1" "${file}"; then
				echo "${file}";
		    	sed "/$2/d;" "${file}" > tmp.dat;
				mv "${file}" "${file}".${timestamp}.bak;
				mv tmp.dat "${file}";
			fi
		done
		;;
	Revert)
		# REVERT CHANGE
		find $1 -name '*.inpt.'"$2.bak" | while read file
		do 
			echo $file
			mv "${file}" "${file%%.${2}.bak}"
		done
		;;
	Test)
		# SUBSTITUTE OR DELETE PATTERN
		find $1 -name '*.inpt' | while read file
		do
		  if grep -q "$2" "${file}"; then
			 echo "${file}";
			 sed s/"$2"/"$3"/g "${file}"
		  fi
		done
		;;
	Search)
		find $1 -name '*.inpt' | while read file
		do
  			if grep "$2" "${file}"; then
				echo "${file}";
			fi
		done
		;;
	*)
		# SUBSTITUTE OR DELETE PATTERN
		find $1 -name '*.inpt' | while read file
		do
		  if grep -q "$2" "${file}"; then
			 echo "${file}";
			 sed s/"$2"/"$3"/g "${file}" > tmp.dat;
			 mv "${file}" "${file}".${timestamp}.bak;
			 mv tmp.dat "${file}";
		  fi
		done
		# Keep history of changes in this file?
		echo "#$1 $2 $3" >> ~/bin/inpt_history.log
		;;
esac

# DELETE SET OF LINES
# find $1 -name '*.inpt' | while read file
# do
#   if grep "$2" "${file}"; then
#      echo "${file}";
#      sed '/$2/$3/}/d' "${file}" > tmp.dat;
#      mv "${file}" "${file}".${timestamp}.bak;
#      mv tmp.dat "${file}";
#   fi
# done
# exit
