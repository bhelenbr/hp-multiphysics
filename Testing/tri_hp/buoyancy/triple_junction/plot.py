#! /usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import math

xle = np.loadtxt("xle.dat", delimiter=' ', unpack=True)
sx = np.loadtxt("sx.dat", delimiter=' ', unpack=True)


plt.plot(sx,xle,'r-x')
plt.xlabel('sx')
plt.ylabel('xle')
plt.savefig("turning.pdf")