#!/usr/bin/env python

import sys
import os
import string
import math
import numpy as np
import matplotlib.pyplot as plt
import glob
from pathlib import Path

os.chdir(os.path.dirname(sys.argv[0]))

# Directory to search
base_dir = "Results"

# load highest order solution and do Richardson Extrapolation
data0 = np.loadtxt(f"{base_dir}/log2p0/cnvg2.dat", delimiter=" ", skiprows=1);
data1 = np.loadtxt(f"{base_dir}/log2p1/cnvg2.dat", delimiter=" ", skiprows=1);
data2 = np.loadtxt(f"{base_dir}/log2p2/cnvg2.dat", delimiter=" ", skiprows=1);

# Just use drag data
data0 = data0[:,[0,1,4]]
data1 = data1[:,[0,1,4]]
data2 = data2[:,[0,1,4]]  

# Exact extrapolated value
alphas = np.log((data2[-2,1:3]-data2[-1,1:3])/(data2[-3,1:3]-data2[-2,1:3]))/math.log(0.5)
cs = (data2[-2,1:3]-data2[-1,1:3])/(2**alphas -1)
exact = data2[-1,1:3]-cs

print(alphas)
print(exact)

errors0 = np.abs(data0[:,1:3] -exact)
errors1 = np.abs(data1[:,1:3] -exact)
errors2 = np.abs(data2[:,1:3] -exact)

resolutions = np.array(range(len(data0)))
resolutions = 2**resolutions
plt.loglog(resolutions,errors0[:,0],'r-x')
plt.loglog(2*resolutions,errors1[:,0],'b-x')
plt.loglog(4*resolutions,errors2[:,0],'g-x')

plt.loglog(resolutions,errors0[:,1],'r--x')
plt.loglog(2*resolutions,errors1[:,1],'b--x')
plt.loglog(4*resolutions,errors2[:,1],'g--x')
		
plt.xlabel('Resolution')
plt.ylabel('Drag Error')
plt.savefig(f"{base_dir}/Error.pdf")
plt.close()	
		
print(np.log2(errors0[-2,:]/errors0[-1,:]))
print(np.log2(errors1[-2,:]/errors1[-1,:]))
print(np.log2(errors2[-2,:]/errors2[-1,:]))


