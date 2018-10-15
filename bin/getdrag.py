#! /usr/bin/env python

import sys
import os
import string
import math
#import numpy
import glob
import fileinput
import re

#os.chdir(os.path.dirname(sys.argv[0]))

filelist = glob.glob("*.drag")

dragfile = open(filelist[0])
i = 0
while 1:
	line = dragfile.readline()
	if not line:
		break
	if (re.search('viscous/pressure flux', line)):
		i+=1

#drag=numpy.zeros((i,2),dtype=float)
drag = [[0.0 for col in range(2)] for row in range(i)]
 
# fill it with 1 diagonally
for i in range(10):
    matrix10x10[i][i] = 1
 
# show it
for row in matrix10x10:
    print row
    
    
for file in filelist:
	dragfile = open(file)
	i = 0
	while 1:
		line = dragfile.readline()
		if not line:
			break
		if (re.search('s7 viscous/pressure flux', line)):
			line = dragfile.readline()
			linelist= string.split(line)
			#drag[i,0] += float(linelist[1])
			#drag[i,1] += float(linelist[2])
			drag[i][0] += float(linelist[1])
			drag[i][1] += float(linelist[2])
			i+=1
	dragfile.close()

for row in drag:
    print row

