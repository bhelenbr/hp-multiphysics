#!/usr/local/bin/python3

import sys
import os
import string
import math
import numpy
import matplotlib.pyplot as plt
import glob

os.chdir(os.path.dirname(sys.argv[0]))

errors = numpy.loadtxt("Results/cnvg.dat", delimiter=" ", skiprows=0);

print(errors)

resolutions = numpy.array([1,2,4,8])

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
plt.savefig("Results/Error.pdf")	

print(math.log2(errors[1,0]/errors[2,0]))
print(math.log2(errors[1,1]/errors[2,1]))
print(math.log2(errors[1,2]/errors[2,2]))
print(math.log2(errors[1,3]/errors[2,3]))