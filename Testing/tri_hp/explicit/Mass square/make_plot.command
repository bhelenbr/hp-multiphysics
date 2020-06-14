#!/usr/local/bin/python3

import sys
import os
#import popen2
import string
import math
import numpy
import matplotlib.pyplot as plt
import glob

os.chdir(os.path.dirname(sys.argv[0]))

errors0 = numpy.loadtxt("Results/cnvg0.dat", delimiter=" ", skiprows=0);
errors1 = numpy.loadtxt("Results/cnvg1.dat", delimiter=" ", skiprows=0);
errors2 = numpy.loadtxt("Results/cnvg2.dat", delimiter=" ", skiprows=0);

resolutions = numpy.array([1,2,4,8,16,32,64])

# L2 errors approximate mass
plt.loglog(resolutions,errors0[0::2,0],'r-x')
plt.loglog(resolutions,errors1[0::2,0],'b-x')
plt.loglog(resolutions,errors2[0::2,0],'g-x')

# L2 errors actual mass
plt.loglog(resolutions,errors0[1::2,0],'r-o')
plt.loglog(resolutions,errors1[1::2,0],'b-o')
plt.loglog(resolutions,errors2[1::2,0],'g-o')

# plt.xlim(rlist[0],rlist[len(rlist)-1])
# plt.ylim(1e-12,1e-2)
# plt.xticks(rlist,rlist)
plt.xlabel('Resolution')
plt.ylabel('Error')
plt.savefig("Results/Error.pdf")	

#print(math.log2(errors0[11,0]/errors0[13,0]))
print(math.log2(errors0[10,0]/errors0[12,0]))
#print(math.log2(errors1[11,0]/errors1[13,0]))
print(math.log2(errors1[10,0]/errors1[12,0]))
#print(math.log2(errors2[11,0]/errors2[13,0]))
print(math.log2(errors2[10,0]/errors2[12,0]))