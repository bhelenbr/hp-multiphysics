#!/usr/bin/env python

import sys
import os
import string
import math
import numpy
import matplotlib.pyplot as plt
import glob

os.chdir(os.path.dirname(sys.argv[0]))

errors0 = numpy.loadtxt("Results/log2p0/cnvg.dat", delimiter=" ", skiprows=0);
errors1 = numpy.loadtxt("Results/log2p1/cnvg.dat", delimiter=" ", skiprows=0);
errors2 = numpy.loadtxt("Results/log2p2/cnvg.dat", delimiter=" ", skiprows=0);
nres = errors0.shape[0]
resolutions = numpy.array(2.0**numpy.array(range(nres)))


# L2 errors
plt.loglog(resolutions,errors0[0::,2],'r-x')
plt.loglog(2*resolutions,errors1[0::,2],'b-x')
plt.loglog(4*resolutions,errors2[0::,2],'g-x')

# Linf errors
plt.loglog(resolutions,errors0[0::,3],'r-o')
plt.loglog(2*resolutions,errors1[0::,3],'b-o')
plt.loglog(4*resolutions,errors2[0::,3],'g-o')

# plt.xlim(rlist[0],rlist[len(rlist)-1])
# plt.ylim(1e-12,1e-2)
# plt.xticks(rlist,rlist)
plt.xlabel('Resolution')
plt.ylabel('Error')
plt.savefig("Results/EvR.pdf")
plt.close()	

# L2 errors vs NDOF
plt.loglog(errors0[0::,0],errors0[0::,2],'r-x')
plt.loglog(errors1[0::,0],errors1[0::,2],'b-x')
plt.loglog(errors2[0::,0],errors2[0::,2],'g-x')

# Linf errors vs NDOF
plt.loglog(errors0[0::,0],errors0[0::,3],'r-o')
plt.loglog(errors1[0::,0],errors1[0::,3],'b-o')
plt.loglog(errors2[0::,0],errors2[0::,3],'g-o')

# plt.xlim(rlist[0],rlist[len(rlist)-1])
# plt.ylim(1e-12,1e-2)
# plt.xticks(rlist,rlist)
plt.xlabel('NDOF')
plt.ylabel('Error')
plt.savefig("Results/EvN.pdf")
plt.close()	

# L2 errors vs cput time
plt.loglog(errors0[0::,1],errors0[0::,2],'r-x')
plt.loglog(errors1[0::,1],errors1[0::,2],'b-x')
plt.loglog(errors2[0::,1],errors2[0::,2],'g-x')

# Linf errors vs NDOF
plt.loglog(errors0[0::,1],errors0[0::,3],'r-o')
plt.loglog(errors1[0::,1],errors1[0::,3],'b-o')
plt.loglog(errors2[0::,1],errors2[0::,3],'g-o')

# plt.xlim(rlist[0],rlist[len(rlist)-1])
# plt.ylim(1e-12,1e-2)
# plt.xticks(rlist,rlist)
plt.xlabel('cpu time')
plt.ylabel('Error')
plt.savefig("Results/EvT.pdf")
plt.close()	

print(math.log2(errors0[nres-2,0]/errors0[nres-1,0]))
print(math.log2(errors1[nres-2,0]/errors1[nres-1,0]))
print(math.log2(errors2[nres-2,0]/errors2[nres-1,0]))

print(math.log2(errors0[nres-2,1]/errors0[nres-1,1]))
print(math.log2(errors1[nres-2,1]/errors1[nres-1,1]))
print(math.log2(errors2[nres-2,1]/errors2[nres-1,1]))