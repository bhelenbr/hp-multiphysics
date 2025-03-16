#! /usr/bin/env python3

import os,sys
import numpy as np
import matplotlib.pyplot as plt
import math

os.chdir(os.path.dirname(sys.argv[0]))

for nrefine in range(1,3):
	cpus, times, jacobian, inversion  = np.loadtxt("Results/nrefine" +str(nrefine) +"/cputimes.dat", delimiter=' ', unpack=True)
	
	nentries = times.size
	theory = np.arange(nentries)+1
	speedup = times[0]/times[:]
	jspeedup = jacobian[0]/jacobian[:]
	ispeedup = inversion[0]/inversion[:]
	
	plt.plot(cpus,speedup,'r-x',label='total' +str(nrefine))
	plt.plot(cpus,jspeedup,'g-x',label='jacobian' +str(nrefine))
	plt.plot(cpus,ispeedup,'b-x',label='petsc inversion' +str(nrefine))
	plt.plot(cpus,theory,'k-',label='theory' +str(nrefine))

plt.legend(loc="upper left")
plt.xlabel('N processor')
plt.ylabel('Speedup')
plt.savefig("Results/speedup.pdf")

