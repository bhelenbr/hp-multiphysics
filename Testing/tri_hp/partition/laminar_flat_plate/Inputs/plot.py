#! /usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt
import math

cpus, times, jacobian, inversion  = np.loadtxt("cputimes.dat", delimiter=' ', unpack=True)

nentries = times.size
theory = np.arange(nentries)+1
speedup = times[0]/times[:]
jspeedup = jacobian[0]/jacobian[:]
ispeedup = inversion[0]/inversion[:]

plt.plot(cpus,speedup,'r-x',label='total')
plt.plot(cpus,jspeedup,'g-x',label='jacobian')
plt.plot(cpus,ispeedup,'b-x',label='petsc inversion')
plt.plot(cpus,theory,'k-',label='theory')
plt.legend(loc="upper left")
plt.xlabel('N processor')
plt.ylabel('Speedup')
plt.savefig("speedup.pdf")
