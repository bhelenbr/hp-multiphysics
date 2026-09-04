#! /usr/bin/env python3

# Simulate an airfoil with a singular point mesh at trailing edge and around airfoil
import sys
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
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

# copy input files into results directory
os.system("cp ../Inputs/* .")

# find intersection points
delta = 0.05
os.system(f"mod_map generate.inpt delta {delta}")

# calculate slope at trailing edge
os.system("spline -i points.dat naca.spl > interp.dat")
s, x, y, tx, ty, curvx, curvy = np.loadtxt("interp.dat", delimiter=' ', unpack=True)

thetaBot = np.atan2(ty[0],tx[0])
thetaTop = np.atan2(ty[1],tx[1])+math.pi
os.system(f"mod_map generate.inpt thetaTop {thetaTop:.8f}")
os.system(f"mod_map generate.inpt thetaBot {thetaBot:.8f}")


# Find top intersection points
s = 2.0-delta

for i in range(100):
	os.system(f"echo {s:.8f} > s.dat")
	os.system("spline -i s.dat naca.spl > interp.dat")
	ss, x, y, tx, ty, curvx, curvy = np.loadtxt("interp.dat", delimiter=' ', unpack=True)
	r2 = (x-1.0)**2+y**2
	dr2ds = 2*(x-1.0)*tx +2*y*ty
	ds = -(r2-delta**2)/dr2ds
	if (math.fabs(ds) < 1.0e-8):
		break
	s = s +ds
	print(s)

sTop = s
xTop = x
yTop = y
theta2 = np.atan2(ty,tx)+math.pi

s = delta
for i in range(100):
	os.system(f"echo {s:.8f} > s.dat")
	os.system("spline -i s.dat naca.spl > interp.dat")
	ss, x, y, tx, ty, curvx, curvy = np.loadtxt("interp.dat", delimiter=' ', unpack=True)
	r2 = (x-1.0)**2+y**2
	dr2ds = 2*(x-1.0)*tx +2*y*ty
	ds = -(r2-delta**2)/dr2ds
	if (math.fabs(ds) < 1.0e-8):
		break
	s = s +ds
	print(s)

sBot = s
xBot = x
yBot = y
theta3 = np.atan2(ty,tx)

os.system(f"mod_map generate.inpt xBot {xBot:.8f}")
os.system(f"mod_map generate.inpt yBot {yBot:.8f}")
os.system(f"mod_map generate.inpt sBot {sBot:.8f}")	
os.system(f"mod_map generate.inpt theta3 {theta3:.8f}")	

os.system(f"mod_map generate.inpt xTop {xTop:.8f}")
os.system(f"mod_map generate.inpt yTop {yTop:.8f}")
os.system(f"mod_map generate.inpt sTop {sTop:.8f}")	
os.system(f"mod_map generate.inpt theta2 {theta2:.8f}")	

os.system("tri_mesh generate.inpt")

os.system("cp generate.inpt run.inpt")
os.system("mod_map run.inpt b0_mesh rstrt3_b0.grd")
os.system("mod_map run.inpt b1_mesh rstrt3_b1.grd")
os.system("mod_map run.inpt \"growth factor\" 10")
os.system("mod_map run.inpt b0_type ins")
os.system("mod_map run.inpt b1_type ins")
os.system("mod_map run.inpt nblock \"1 1\"")
os.system("mod_map run.inpt logfile output")
os.system("mod_map run.inpt adapt 0")
os.system("mod_map run.inpt ncycle 10")
os.system("mod_map run.inpt ntstep 1")
#os.system("mod_map run.inpt rsdl_debug 1")

os.system("mpiexec -np 2 tri_hp_petsc run.inpt")

os.system("mod_map run.inpt log2p 1")
os.system("mod_map run.inpt restart 1")
os.system("mpiexec -np 2 tri_hp_petsc run.inpt")

os.system("mod_map run.inpt log2p 2")
os.system("mod_map run.inpt restart 2")
os.system("mpiexec -np 2 tri_hp_petsc run.inpt")



os.chdir("..")
os.system("opendiff Baseline/ Results/")
#os.system("open naca.pdf");