#!/usr/bin/env python

import sys
import os
import string
import math
import numpy
import matplotlib.pyplot as plt
import glob

os.chdir(os.path.dirname(sys.argv[0]))

filename = "Results/data0_b0_s1.dat"
freesurface = numpy.loadtxt(filename, skiprows=3)
plt.plot(freesurface[:,1],freesurface[:,2],'b-x')

filename = "Results/data1_b0_s1.dat"
freesurface = numpy.loadtxt(filename, skiprows=3)
plt.plot(freesurface[:,1],freesurface[:,2],'b-x')	

filename = "Results/data2_b0_s1.dat"
freesurface = numpy.loadtxt(filename, skiprows=3)
plt.plot(freesurface[:,1],freesurface[:,2],'b-x')

filename = "Results/data3_b0_s7.dat"
freesurface = numpy.loadtxt(filename, skiprows=3)
plt.plot(freesurface[:,1],freesurface[:,2],'g-x')

filename = "Results/data3_b2_s12.dat"
freesurface = numpy.loadtxt(filename, skiprows=3)
plt.plot(freesurface[:,1],freesurface[:,2],'r-x')

filename = "Results/data4_b0_s7.dat"
freesurface = numpy.loadtxt(filename, skiprows=3)
plt.plot(freesurface[:,1],freesurface[:,2],'g-x')

filename = "Results/data4_b2_s12.dat"
freesurface = numpy.loadtxt(filename, skiprows=3)
plt.plot(freesurface[:,1],freesurface[:,2],'r-x')

plt.xlabel('x')
plt.ylabel('y')
plt.savefig("Results/free_surface.pdf")	