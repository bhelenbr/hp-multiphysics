#!/usr/bin/env python

import sys
import os
import string
import math
import numpy
import matplotlib.pyplot as plt
import glob

os.chdir(os.path.dirname(sys.argv[0]))

errors0 = numpy.loadtxt("Results_petsc/cnvg.dat", delimiter=" ", skiprows=0);

n = errors0.shape[0]
resolutions = numpy.arange(0,n)
resolutions = 2**resolutions

plt.loglog(resolutions,errors0[0::,0],'r-x')
plt.loglog(resolutions,errors0[0::,1],'r--x')
plt.xlabel('Resolution')
plt.ylabel('Error')
plt.savefig("Results_petsc/Error.pdf")
plt.close()	

print('L2 convergence')
print(math.log2(errors0[n-2,0]/errors0[n-1,0]))
print(math.log2(errors0[n-2,1]/errors0[n-1,1]))
