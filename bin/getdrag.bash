#!/bin/bash

# If the files have different numbers of entries this won't work
# They should all have the same, but maybe there should be a test using 'wc -l'


# For each log file containing the word "viscous"
# Make a new file with x,y forces in columns
for file in $(grep -l viscous *.log); do
	echo $file
	grep -A 1 viscous $file > ${file%.log}.drag;
done

# Process the drag files using python
getdrag.py
