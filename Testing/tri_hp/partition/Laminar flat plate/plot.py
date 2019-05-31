#! /usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import math

cpus, times = np.loadtxt("cputimes.dat", delimiter=' ', unpack=True)

nentries = times.size
theory = np.arange(nentries)+1
speedup = times[0]/times[:]

plt.plot(cpus,speedup,'r-x')
plt.plot(cpus,theory,'b-')
plt.xlabel('N processor')
plt.ylabel('Speedup')
plt.savefig("speedup.pdf")
