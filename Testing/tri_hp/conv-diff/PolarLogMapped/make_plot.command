#!/usr/bin/env python

import sys
import os
import string
import math
import numpy
import matplotlib.pyplot as plt
import glob
from pathlib import Path

os.chdir(os.path.dirname(sys.argv[0]))

# Directory to search
base_dir = Path("Results")

# Loop through folders
for folder in base_dir.iterdir():
	if folder.is_dir():
		print(f"Folder: {folder}")
		
		
		errors0 = numpy.loadtxt(f"{folder}/log2p0/cnvg.dat", delimiter=" ", skiprows=0);
		errors1 = numpy.loadtxt(f"{folder}/log2p1/cnvg.dat", delimiter=" ", skiprows=0);
		errors2 = numpy.loadtxt(f"{folder}/log2p2/cnvg.dat", delimiter=" ", skiprows=0);
		resolutions = numpy.array([1,2,4,8])
		
		# L2 errors
		plt.loglog(resolutions,errors0[0::,1],'r-x')
		plt.loglog(2*resolutions,errors1[0::,1],'b-x')
		plt.loglog(4*resolutions,errors2[0::,1],'g-x')
		
		# Linf errors
		plt.loglog(resolutions,errors0[0::,2],'r-o')
		plt.loglog(2*resolutions,errors1[0::,2],'b-o')
		plt.loglog(4*resolutions,errors2[0::,2],'g-o')
		
		# plt.xlim(rlist[0],rlist[len(rlist)-1])
		# plt.ylim(1e-12,1e-2)
		# plt.xticks(rlist,rlist)
		plt.xlabel('Resolution')
		plt.ylabel('Error')
		plt.savefig(f"{folder}/Error.pdf")
		plt.close()	
		
		print(math.log2(errors0[2,1]/errors0[3,1]))
		print(math.log2(errors1[2,1]/errors1[3,1]))
		print(math.log2(errors2[2,1]/errors2[3,1]))
		
		print(math.log2(errors0[2,2]/errors0[3,2]))
		print(math.log2(errors1[2,2]/errors1[3,2]))
		print(math.log2(errors2[2,2]/errors2[3,2]))\
		
		# L2 errors
		plt.loglog(errors0[0::,0],errors0[0::,1],'r-x')
		plt.loglog(errors1[0::,0],errors1[0::,1],'b-x')
		plt.loglog(errors2[0::,0],errors2[0::,1],'g-x')
		
		# Linf errors
		plt.loglog(errors0[0::,0],errors0[0::,2],'r-o')
		plt.loglog(errors1[0::,0],errors1[0::,2],'b-o')
		plt.loglog(errors2[0::,0],errors2[0::,2],'g-o')
		
		# plt.xlim(rlist[0],rlist[len(rlist)-1])
		# plt.ylim(1e-12,1e-2)
		# plt.xticks(rlist,rlist)
		plt.xlabel('DOF')
		plt.ylabel('Error')
		plt.savefig(f"{folder}/ErrorDOF.pdf")
		plt.close()	


