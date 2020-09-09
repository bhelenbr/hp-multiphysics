#! /usr/bin/env python

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

plt.plot(cpus,speedup,'r-x')
plt.plot(cpus,jspeedup,'g-x')
plt.plot(cpus,ispeedup,'b-x')
plt.plot(cpus,theory,'k-')
plt.xlabel('N processor')
plt.ylabel('Speedup')
plt.savefig("speedup.pdf")
