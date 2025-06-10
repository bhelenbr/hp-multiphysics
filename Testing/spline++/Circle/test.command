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
os.environ['PATH'] = BINDIR.strip().decode('ascii') + os.environ['PATH']

if not os.path.isdir("Results"):
	os.mkdir("Results")
os.chdir("Results")
os.system("rm *")

error = []
nsegs = 4
nres = 5
xy = plt.figure()
tan = plt.figure()
crv = plt.figure()
for cnt in range(nres):
	with open("circle.spl", "w") as file:
		file.write("Circle Surface (SCALED 0-2pi)\n")
		file.write("NPTS: " +str(nsegs+1) +"\n")
		file.write("S-COORD          X                Y\n")
		for i in range(nsegs+1):
			angle = 2 * math.pi * i / nsegs
			x = math.cos(angle)
			y = math.sin(angle)
			file.write(f"{angle:0.16f} {x:.16f} {y:.16f}\n")
	
	theta_vals = []
	with open("pts.dat", "w") as file:
		nsegs=4*nsegs
		for i in range(nsegs+1):
			angle = 2 * math.pi * i / nsegs
			file.write(f"{angle:0.16f}\n")
			theta_vals.append(angle)
	
	os.system("spline -i pts.dat circle.spl > circle.dat")
	
	s_vals = []
	x_vals = []
	y_vals = []
	xtan_vals = []
	ytan_vals = []
	xcrv_vals = []
	ycrv_vals = []
	
	with open('circle.dat', 'r') as file:
		for line in file:
			if line.strip():  # skip empty lines
				parts = line.strip().split()
				if len(parts) >= 6:
					s_vals.append(float(parts[0]))  # 2nd column
					x_vals.append(float(parts[1]))  # 2nd column
					y_vals.append(float(parts[2]))  # 3rd column
					xtan_vals.append(float(parts[3]))  # 2nd column
					ytan_vals.append(float(parts[4]))  # 2nd column
					xcrv_vals.append(float(parts[5]))  # 2nd column
					ycrv_vals.append(float(parts[6]))  # 2nd column

					
	x = numpy.array(x_vals)
	y = numpy.array(y_vals)
	theta = numpy.array(theta_vals)
	radii = numpy.abs(numpy.sqrt(x**2 + y**2)-1)
	
	error.append(numpy.sum(radii))
	nsegs = nsegs*2

	# Plotting
	plt.figure(xy)
	plt.plot(x_vals, y_vals, 'o-')
	
	plt.figure(tan)
	plt.plot(s_vals,xtan_vals,'r')
	plt.plot(s_vals,ytan_vals,'b')
	plt.plot(theta_vals,-numpy.sin(theta),'r--')
	plt.plot(theta_vals,numpy.cos(theta),'b--')
	
	plt.figure(crv)
	plt.plot(s_vals,xcrv_vals,'r')
	plt.plot(s_vals,ycrv_vals,'b')
	plt.plot(theta_vals,-numpy.cos(theta),'r--')
	plt.plot(theta_vals,-numpy.sin(theta),'b--')	
	

plt.figure(xy)	
plt.grid(True)
plt.savefig('circle.pdf')
plt.close()

plt.figure(tan)	
plt.grid(True)
plt.savefig('tan.pdf')
plt.close()

plt.figure(crv)	
plt.grid(True)
plt.savefig('crv.pdf')
plt.close()


resolutions = numpy.array(4*2.0**numpy.array(range(nres)))

# L2 errors
plt.figure()
plt.loglog(resolutions,error,'r-x')
plt.savefig('error.pdf')

with open('error.dat', 'w') as file:
	for i in range(nres):
		file.write(f"{resolutions[i]:0.16f} {error[i]:0.16f}\n")


os.chdir('..')
os.system('opendiff Results/ Baseline/')
