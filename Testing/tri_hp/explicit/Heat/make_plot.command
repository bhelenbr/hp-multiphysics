#!/usr/bin/python3

import sys
import os
import string
import math
import numpy
import matplotlib.pyplot as plt
import glob

os.chdir(os.path.dirname(sys.argv[0]))

errors0 = numpy.loadtxt("Results/cnvg0.dat", delimiter=" ", skiprows=0);
errors1 = numpy.loadtxt("Results/cnvg1.dat", delimiter=" ", skiprows=0);
errors2 = numpy.loadtxt("Results/cnvg2.dat", delimiter=" ", skiprows=0);
errors3 = numpy.loadtxt("Results/cnvg3.dat", delimiter=" ", skiprows=0);

resolutions = numpy.array([1,2,4])

# L2 errors approximate mass
plt.loglog(resolutions,errors0[0::,0],'r-x')
plt.loglog(resolutions,errors1[0::,0],'b-x')
plt.loglog(resolutions,errors2[0::,0],'g-x')
plt.loglog(resolutions,errors3[0::,0],'k-x')

# plt.xlim(rlist[0],rlist[len(rlist)-1])
# plt.ylim(1e-12,1e-2)
# plt.xticks(rlist,rlist)
plt.xlabel('Resolution')
plt.ylabel('Error')
plt.savefig("Results/Error.pdf")	

print(math.log2(errors0[1,0]/errors0[2,0]))
print(math.log2(errors1[1,0]/errors1[2,0]))
print(math.log2(errors2[1,0]/errors2[2,0]))
print(math.log2(errors3[1,0]/errors3[2,0]))