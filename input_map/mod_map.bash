#!/bin/bash

while getopts ":dcuea" opt; do
   case $opt in
      d) Command="Delete";;
      c) Command="Comment";;
      u) Command="Uncomment";;
      e) Command="Echo";;
      a) Command="Append";;
      \?) echo 'usage\: mod_map.script -[d,c,u]'
           exit 1
   esac
done
shift $((OPTIND-1))


# Using :'s as delimiter for replace string
# because it is already a special character anyway
# avoids problems with replace strings with a /
case "${Command}" in 
	Delete)
		sed -i.bak -E "s:^ *$2 *\:.*$::g" $1;;
	Comment)
		sed -i.bak -E "s:^ *$2 *\::#$2\::g" $1;;
	Uncomment)
		sed -i.bak -E "s:^# *$2 *\::$2\::g" $1;;
	Echo)
		grep "^ *$2 *:.*$" $1 | tail -n 1 | cut -f 2 -d: | sed "s:^ *::";;
	*)
		if grep -q "^ *$2 *:.*$" $1; then
			sed -i.bak -E "s:^ *$2 *\:.*$:$2\: $3:g" $1;
		else
			echo "$2: $3" >> $1
		fi
esac