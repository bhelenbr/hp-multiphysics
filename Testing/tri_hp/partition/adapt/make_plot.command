#!/usr/bin/env python

import sys
import os
import string
import math
import numpy
import matplotlib.pyplot as plt
import glob

os.chdir(os.path.dirname(sys.argv[0]))

filename = "Results/data0_b0_s3.dat"
freesurface = numpy.loadtxt(filename, skiprows=3)
plt.plot(freesurface[:,1],freesurface[:,2],'b-x')

filename = "Results/data1_b0_s3.dat"
freesurface = numpy.loadtxt(filename, skiprows=3)
plt.plot(freesurface[:,1],freesurface[:,2],'g-x')	
	
for block in range(3):
	filename = "Results/data2_b" +str(block) +"_s3.dat"
	freesurface = numpy.loadtxt(filename, skiprows=3)
	plt.plot(freesurface[:,1],freesurface[:,2],'r-x')

plt.xlabel('x')
plt.ylabel('y')
plt.savefig("Results/free_surface.pdf")	