#!/bin/sh
##
# script by Brian Helenbrook | helenbrk@clarkson.edu
##

location=$(dirname "$0")
for file; do
    case $file in
      b?_*.dat)  
			cat $file >> ${file#b?_};
			rm $file;;
    esac
done
