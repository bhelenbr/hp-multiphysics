#!/bin/bash

timestamp=`date +%Y%m%d%H%M`

# find . -name '*.command' | while read file
# do
#   if grep -q "delete_data.bash" "${file}"; then
# 	 echo "${file}";
# 	 sed 's/delete_data.bash \([0-9]*\) \([0-9]*\)/delete_series data \1 \2\ \
# > delete_series rstrt 1 1/g' "${file}" > tmp.dat;
# 	 mv "${file}" "${file}".${timestamp}.bak;
# 	 mv tmp.dat "${file}";
#   fi
# done

find . -name '*.command' | while read file
do
  if grep -q 'delete_series' "${file}"; then
	 echo "${file}";
	 sed 's/delete_series data 2 0/delete_series data 0 1/g' "${file}" > tmp.dat;
	 mv "${file}" "${file}".${timestamp}.bak;
	 mv tmp.dat "${file}";
	 chmod +x "${file}"
  fi
done
