#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os
import sys

os.chdir(os.path.dirname(sys.argv[0]))
data91 = np.loadtxt('Results91/turning.dat',delimiter=' ')
data10 = np.loadtxt('Results10/turning.dat',delimiter=' ')
data5 = np.loadtxt('Results5/turning.dat',delimiter=' ')

first_zero = np.argmax(data91[:,0]==0)
data91 = data91[:first_zero-1,:]
first_zero = np.argmax(data10[:,0]==0)
data10 = data10[:first_zero-1,:]
first_zero = np.argmax(data5[:,0]==0)
data5 = data5[:first_zero-1,:]

plt.plot(data91[:,0],data91[:,1])
plt.plot(data10[:,0],data10[:,1])
plt.plot(data5[:,0],data5[:,1])

plt.xlabel('sx')
plt.ylabel('xle')
plt.savefig('compare.pdf')
plt.close()	 

