#!/usr/bin/env python

import sys
import os
import string
import math
import numpy
import matplotlib.pyplot as plt
import glob

#os.chdir(os.path.dirname(sys.argv[0]))
errors = numpy.loadtxt(sys.argv[1] +"/cnvg.dat", delimiter=" ", skiprows=0);
n = errors.shape[0]

#print(errors)

resolutions = 2**numpy.arange(0,n)

# L2 errors approximate mass
plt.loglog(resolutions,errors[0::,0],'r-x')
plt.loglog(resolutions,errors[0::,1],'b-x')
plt.loglog(resolutions,errors[0::,2],'g-x')
plt.loglog(resolutions,errors[0::,3],'k-x')

# plt.xlim(rlist[0],rlist[len(rlist)-1])
# plt.ylim(1e-12,1e-2)
# plt.xticks(rlist,rlist)
plt.xlabel('Resolution')
plt.ylabel('Error')
plt.savefig(sys.argv[1] +"/Error.pdf")	

print(math.log2(errors[n-2,0]/errors[n-1,0]))
print(math.log2(errors[n-2,1]/errors[n-1,1]))
print(math.log2(errors[n-2,2]/errors[n-1,2]))
print(math.log2(errors[n-2,3]/errors[n-1,3]))