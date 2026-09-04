#!/usr/bin/env python3

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


def run(cmd, cwd=None):
    """Run a command and raise an exception on failure."""
    print("Running:", " ".join(map(str, cmd)))
    subprocess.run(cmd, cwd=cwd, check=True)

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

offset = -0.05
os.system("spline -m 0.0,0.0 -r 0.0 -s 1.0 -o" +str(offset)+ " -i spoints.dat naca.spl > interp.dat");
s, x, y, tx, ty, curvx, curvy = np.loadtxt("interp.dat", delimiter=' ', unpack=True)
npoints = len(s)

for n in range(5):
	os.system(f"mod_map generate.inpt x{n} {x[n]:.16f}")
	os.system(f"mod_map generate.inpt y{n} {y[n]:.16f}")
	os.system(f"mod_map generate.inpt s{n} {s[n]:.16f}")
os.system("tri_mesh generate.inpt")

os.system("cp generate.inpt run.inpt")
os.system("mod_map run.inpt b0_mesh rstrt0_b0.grd")
os.system("mod_map run.inpt b1_mesh rstrt0_b1.grd")
os.system("mod_map run.inpt b2_mesh rstrt0_b2.grd")
os.system("mod_map run.inpt \"growth factor\" 10")
os.system("mod_map run.inpt b0_type ins")
os.system("mod_map run.inpt b1_type ins")
os.system("mod_map run.inpt b2_type ins")
os.system("mod_map run.inpt nblock \"1 1 1\"")
os.system("mod_map run.inpt logfile output")
os.system("mod_map run.inpt adapt 0")
os.system("mod_map run.inpt ncycle 10")
os.system("mod_map run.inpt ntstep 1")

for log2p in range(1,3):
	case_dir = Path(f"log2p{log2p}")
	case_dir.mkdir()

	shutil.copy("run.inpt", case_dir)
	shutil.copy("naca.spl", case_dir)
	shutil.copy("rstrt3_b0.grd", case_dir / "rstrt0_b0.grd")
	shutil.copy("rstrt3_b1.grd", case_dir / "rstrt0_b1.grd")
	shutil.copy("rstrt3_b2.grd", case_dir / "rstrt0_b2.grd")
	os.chdir(case_dir)

	run(["mod_map", "run.inpt", "log2p", str(log2p)])

	run(["mpiexec", "-np", "3", "tri_hp_petsc", "run.inpt"])
	
	ngrids = 3
	restart = 1

	for ngrid in range(1, ngrids):
	
		# Refine solution
		run(["mod_map", "run.inpt", "refineby2", "1"])
		run(["mod_map", "run.inpt", "adapt", "1"])
		run(["mod_map", "run.inpt", "restart", str(restart)])
	
		run(["mpiexec", "-np", "3", "tri_hp_petsc", "run.inpt"])
		restart += 1
	
		# Run restart case
		run(["mod_map", "run.inpt", "restart", str(restart)])
		run(["mod_map", "run.inpt", "adapt", "0"])
		run(["mod_map", "run.inpt", "refineby2", "0"])
		run(["mod_map", "run.inpt", "b0_mesh", f"rstrt{restart}_b0.nc"])
		run(["mod_map", "run.inpt", "b1_mesh", f"rstrt{restart}_b1.nc"])
		run(["mod_map", "run.inpt", "b2_mesh", f"rstrt{restart}_b2.nc"])
		run(["mpiexec", "-np", "3", "tri_hp_petsc", "run.inpt"])
		restart += 1	
	os.chdir("..")

os.chdir("..")
run(["opendiff", "Results/", "Baseline/"])
