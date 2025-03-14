#! /usr/bin/env python

import sys
import os
import string
import math
import numpy
import matplotlib.pyplot as plt
import glob

os.chdir(os.path.dirname(sys.argv[0]))

dragdns = numpy.loadtxt("Results/DNS/drag.dat", delimiter=" ", skiprows=1);
dragpod = numpy.loadtxt("Results/POD_SIM/drag.dat", delimiter=" ");
plt.plot(dragdns[:,0],dragdns[:,2],'r-x')
plt.plot(dragpod[:,0],dragpod[:,2],'b-x')

# plt.xlim(rlist[0],rlist[len(rlist)-1])
# plt.ylim(1e-12,1e-2)
# plt.xticks(rlist,rlist)
plt.xlabel('Time Step')
plt.ylabel('Drag')
plt.savefig("Results/Drag.pdf")	