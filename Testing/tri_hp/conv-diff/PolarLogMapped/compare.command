#!/usr/bin/env python

import sys
import os
import string
import math
import numpy
import matplotlib.pyplot as plt
import glob

os.chdir(os.path.dirname(sys.argv[0]))

#errors0 = numpy.loadtxt("PolarLogMapped/Results/log2p0/cnvg.dat", delimiter=" ", skiprows=0);
merrors1 = numpy.loadtxt("Results/log2p1/cnvg.dat", delimiter=" ", skiprows=0);
merrors2 = numpy.loadtxt("Results/log2p2/cnvg.dat", delimiter=" ", skiprows=0);

aerrors1 = numpy.loadtxt("../PolarAdaptive/Results/log2p1/cnvg.dat", delimiter=" ", skiprows=0);
aerrors2 = numpy.loadtxt("../PolarAdaptive/Results/log2p2/cnvg.dat", delimiter=" ", skiprows=0);


# L2 errors
#plt.loglog(errors0[0::,2],errors0[0::,0],'r-x')
plt.loglog(aerrors1[0::,2],aerrors1[0::,0],'b--x')
plt.loglog(aerrors2[0::,2],aerrors2[0::,0],'g--x')
plt.loglog(merrors1[0::,2],merrors1[0::,0],'b-x')
plt.loglog(merrors2[0::,2],merrors2[0::,0],'g-x')

# Linf errors
#plt.loglog(errors0[0::,2],errors0[0::,1],'r-o')
plt.loglog(aerrors1[0::,2],aerrors1[0::,1],'b--o')
plt.loglog(aerrors2[0::,2],aerrors2[0::,1],'g--o')
plt.loglog(merrors1[0::,2],merrors1[0::,1],'b-o')
plt.loglog(merrors2[0::,2],merrors2[0::,1],'g-o')

# plt.xlim(rlist[0],rlist[len(rlist)-1])
# plt.ylim(1e-12,1e-2)
# plt.xticks(rlist,rlist)
plt.xlabel('DOF')
plt.ylabel('Error')
plt.savefig("Results/Compare.pdf")
plt.close()	

