#! /usr/bin/env python3

# Uses a file of spline points to create boundary layer mesh points around an airfoil
import sys
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import string
import math

os.chdir(os.path.dirname(sys.argv[0]))

# Define location of executables
p0 = subprocess.Popen("echo ${PWD%/*/*/*}/bin/:", stdout=subprocess.PIPE,shell=True,text=True)
(BINDIR, err) = p0.communicate()
os.environ['PATH'] = BINDIR[:-1] + os.environ['PATH']

if not os.path.isdir("Results"):
	os.mkdir("Results")
os.chdir("Results")
os.system("rm *")

# copy input files into results directory
os.system("cp ../Inputs/* .")
os.system("cp ../test.command .")

# generate mesh and remove unnecessary data files
npoints = 10
offset = 0.0
#os.system("spline -m 0.25,0.0 -r -10.0 -s 3.0 -o" +str(offset)+ " -i spoints.dat naca.spl > interp.dat");
os.system("spline -m 0.0,0.0 -r -0.0 -s 1.0 -o " +str(offset)+ " naca.spl " +str(npoints) + " 3.0 >> interp.dat");
s, x, y, tx, ty, curvx, curvy = np.loadtxt("interp.dat", delimiter=' ', unpack=True)

# plt.plot(x,y,'r-x')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.axis('equal')
# plt.savefig("naca.pdf")

npoints = len(s)
# output .d files
f = open('blayer.d','w')
f.write(str(npoints-1 +4)+'\n')

count = 0
# airfoil trailing edge
f.write('{0:3d}: {1:.10f} {2:.10f} 0.025 1\n'.format(count,x[0],y[0]))
#skip trailing edge point
for i in range(1,npoints-1):
	f.write('{0:3d}: {1:.10f} {2:.10f} 0.025 0\n'.format(count,x[i],y[i]))
	count += 1

#write outer domain points
extent = 20.0
f.write('{0:3d}: {1:.10f} {2:.10f} 2.0 0\n'.format(count,-extent,-extent))
count += 1
f.write('{0:3d}: {1:.10f} {2:.10f} 2.0 0\n'.format(count,extent,-extent))
count += 1
f.write('{0:3d}: {1:.10f} {2:.10f} 2.0 0\n'.format(count,extent,extent))
count += 1
f.write('{0:3d}: {1:.10f} {2:.10f} 2.0 0\n'.format(count,-extent,extent))
count += 1
	
f.write('{0:d}\n'.format(npoints-1 +4))
# Outer boundary layer surface
count = 0
for i in range(npoints-2):
	f.write('{0:d}: {1:d} {2:d} 1\n'.format(count,i,i+1))
	count += 1
# connection to wake point bot
f.write('{0:d}: {1:d} {2:d} 1\n'.format(count,npoints-2,0))
count += 1

# Outer boundary
f.write('{0:d}: {1:d} {2:d} 2\n'.format(count,npoints-1,npoints))
count += 1
f.write('{0:d}: {1:d} {2:d} 3\n'.format(count,npoints,npoints+1))
count += 1
f.write('{0:d}: {1:d} {2:d} 4\n'.format(count,npoints+1,npoints+2))
count += 1
f.write('{0:d}: {1:d} {2:d} 5\n'.format(count,npoints+2,npoints-1))
count += 1

f.close()

os.system("tri_mesh generate.inpt")

os.system("tri_mesh offset.inpt")

os.chdir("..")
#os.system("opendiff Baseline/ Results/")
#os.system("open naca.pdf");