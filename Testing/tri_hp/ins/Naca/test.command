#! /usr/bin/env python3

# Simulate an airfoil with a singular point mesh at trailing edge and around airfoil
import sys
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import string
import math
import shutil
from pathlib import Path

FIELD_SEP = ' '


def cut_fields(line, fields):
	"""Mimic `cut -d' ' -f<fields>` (1-indexed field numbers)."""
	parts = line.rstrip('\n').split(FIELD_SEP)
	selected = [parts[i - 1] for i in fields]
	return FIELD_SEP.join(selected)


def grep_first(filename, needle, offset=0):
	with open(filename) as f:
		lines = f.readlines()
 
	for i, line in enumerate(lines):
		if needle in line:
			target = i + offset
			if target < 0 or target >= len(lines):
				raise ValueError(
					f"offset {offset} from line {i + 1} in {filename} "
					f"falls outside the file (which has {len(lines)} lines)"
				)
			return lines[target]
 
	raise ValueError(f"no line containing {needle!r} found in {filename}")

def grep_last(filename, needle, offset=0):
	with open(filename) as f:
		lines = f.readlines()
 
	for i in range(len(lines) - 1, -1, -1):
		if needle in lines[i]:
			target = i + offset
			if target < 0 or target >= len(lines):
				raise ValueError(
					f"offset {offset} from line {i + 1} in {filename} "
					f"falls outside the file (which has {len(lines)} lines)"
				)
			return lines[target]
 
	raise ValueError(f"no line containing {needle!r} found in {filename}")

def tail(lines, n):
	return lines[-n:]

def run(cmd, cwd=None):
	"""Run a command and raise an exception on failure."""
	print("Running:", " ".join(map(str, cmd)))
	subprocess.run(cmd, cwd=cwd, check=True)
	
def getdata():
	dof = int(cut_fields(grep_last('output_b0.log', "DOF:"), [6]))

	viscous1 = cut_fields(grep_last('output_b0.log','viscous/pressure',1), [2,3,4])
	total1 = cut_fields(grep_last('output_b0.log','total fluxes',1), [2,3,4])
	
	with open('cnvg.dat', "a") as out:
		out.write(f"{dof} {viscous1} {total1}\n")


os.chdir(os.path.dirname(sys.argv[0]))

# Define location of executables
p0 = subprocess.Popen("echo ${PWD%/Testing/*}/bin/:", stdout=subprocess.PIPE,shell=True)
(BINDIR, err) = p0.communicate()
os.environ['PATH'] = BINDIR.strip().decode('ascii') + os.environ['PATH']

if not os.path.isdir("Results"):
	os.mkdir("Results")
os.chdir("Results")
os.system("rm -rf *")

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
os.system("mod_map run.inpt \"growth factor\" 10")
os.system("mod_map run.inpt b0_type ins")
os.system("mod_map run.inpt nblock 1")
os.system("mod_map run.inpt logfile output")
os.system("mod_map run.inpt adapt 0")
os.system("mod_map run.inpt ncycle 10")
os.system("mod_map run.inpt ntstep 1")

os.system("mod_map run.inpt b0_mesh rstrt0_b0.grd")

for log2p in range(3):
	case_dir = Path(f"log2p{log2p}")
	case_dir.mkdir()

	shutil.copy("run.inpt", case_dir)
	shutil.copy("naca.spl", case_dir)
	shutil.copy("rstrt3_b0.grd", case_dir / "rstrt0_b0.grd")
	os.chdir(case_dir)

	run(["mod_map", "run.inpt", "log2p", str(log2p)])

	run(["mpiexec", "-np", "1", "tri_hp_petsc", "run.inpt"])
	getdata()

	ngrids = 4
	restart = 1

	for ngrid in range(1, ngrids):
	
		# Refine solution
		run(["mod_map", "run.inpt", "refineby2", "1"])
		run(["mod_map", "run.inpt", "adapt", "1"])
		run(["mod_map", "run.inpt", "restart", str(restart)])
	
		run(["mpiexec", "-np", "1", "tri_hp_petsc", "run.inpt"])
		restart += 1
	
		# Run restart case
		run(["mod_map", "run.inpt", "restart", str(restart)])
		run(["mod_map", "run.inpt", "adapt", "0"])
		run(["mod_map", "run.inpt", "refineby2", "0"])
		run(["mod_map", "run.inpt", "b0_mesh", f"rstrt{restart}_b0.nc"])
		run(["mpiexec", "-np", "1", "tri_hp_petsc", "run.inpt"])
		getdata()
		restart += 1
	os.chdir("..")
os.chdir("..")
subprocess.run(["python3", "make_plot.command"])
