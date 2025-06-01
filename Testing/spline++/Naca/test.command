#!/usr/bin/env python
# This is a test of spline interpolation for a circle
import sys
import os
import numpy
import matplotlib.pyplot as plt
import subprocess
#import sr
import glob
import string
import math

os.chdir(os.path.dirname(sys.argv[0]))

# Define location of executables
p0 = subprocess.Popen("echo ${PWD%/Testing/*}/bin/:", stdout=subprocess.PIPE,shell=True)
(BINDIR, err) = p0.communicate()
os.environ['PATH'] = str(BINDIR[:-1]) + os.environ['PATH']

if not os.path.isdir("Results"):
	os.mkdir("Results")
os.chdir("Results")
os.system("rm *")

# copy input files into results directory
os.system("cp ../Inputs/* .")

error = []
nsegs = 4000
nres = 5
smax = 0.9
smin = 1.1
omin = -0.07
omax = -0.01

with open("pts.dat", "w") as file:
	for i in range(nsegs+1):
		s = smin +(smax -smin)* i / nsegs
		file.write(f"{s:0.16f}\n")
	
	
noffset = 4
f1 = plt.figure()
f2 = plt.figure()
f3 = plt.figure()
for k in range(noffset+1):
	offset = omin +(omax-omin)*k/noffset
	
	os.system(f"spline -o  {offset:0.10f} -i pts.dat naca.spl > spline.dat")
		
	s_vals = []
	x_vals = []
	y_vals = []
	
	with open('spline.dat', 'r') as file:
		for line in file:
			if line.strip():  # skip empty lines
				parts = line.strip().split()
				if len(parts) >= 3:
					s_vals.append(float(parts[0]))  # 1st column
					x_vals.append(float(parts[1]))  # 2nd column
					y_vals.append(float(parts[2]))  # 3rd column
					
	# Plotting
	plt.figure(f1)
	plt.plot(s_vals, x_vals, '-',label=f"{offset:0.2f}")
	
	plt.figure(f2)
	plt.plot(s_vals, y_vals, '-',label=f"{offset:0.2f}")
	
	plt.figure(f3)
	plt.plot(x_vals, y_vals, '-',label=f"{offset:0.2f}")
	

#Error in hp_edge_bdry::calc_metrics b0_s10 1 (-4.839e-02,-1.443e-02) (1.040e+00,1.185e-02) (-4.971e-02,-6.147e-03) (1.018e+00,5.000e-02)


xypt = [-4.839e-02,-1.443e-02]
sopt = [1.040e+00,1.185e-02]
plt.figure(f1)
plt.plot(sopt[0],xypt[0],'x')
plt.figure(f2)
plt.plot(sopt[0],xypt[1],'x')
plt.figure(f3)
plt.plot(xypt[0],xypt[1],'x')

xypt = [-4.971e-02,-6.147e-03]
sopt = [1.018e+00,5.000e-02]
plt.figure(f1)
plt.plot(sopt[0],xypt[0],'x')
plt.figure(f2)
plt.plot(sopt[0],xypt[1],'x')
plt.figure(f3)
plt.plot(xypt[0],xypt[1],'x')



plt.figure(f1)
plt.grid(True)
plt.legend()
plt.savefig('nacax.pdf')
plt.close()
	
plt.figure(f2)
plt.grid(True)
plt.legend()
plt.savefig('nacay.pdf')
plt.close()

plt.figure(f3)
plt.grid(True)
plt.legend()
plt.savefig('naca.pdf')
plt.close()

os.chdir('..')
os.system('opendiff Results/ Baseline/')
