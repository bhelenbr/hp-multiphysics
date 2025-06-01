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
p0 = subprocess.Popen("echo ${PWD%/Testing/*}/bin/:", stdout=subprocess.PIPE,shell=True)
(BINDIR, err) = p0.communicate()
os.environ['PATH'] = BINDIR[:-1] + os.environ['PATH']

if not os.path.isdir("Results"):
	os.mkdir("Results")
os.chdir("Results")
os.system("rm *")

# copy input files into results directory
os.system("cp ../Inputs/* .")

os.system("tri_mesh generate.inpt")

os.system("cp generate.inpt run.inpt")
os.system("mod_map run.inpt b0_mesh rstrt3_b0.grd")
os.system("mod_map run.inpt b1_mesh rstrt3_b1.grd")
os.system("mod_map run.inpt b2_mesh rstrt3_b2.grd")
os.system("mod_map run.inpt \"growth factor\" 10")
os.system("mod_map run.inpt b0_type ins")
os.system("mod_map run.inpt b1_type ins")
os.system("mod_map run.inpt b2_type ins")
os.system("mod_map run.inpt nblock \"1 1 1\"")
os.system("mod_map run.inpt logfile output")
os.system("mod_map run.inpt adapt 0")
os.system("mod_map run.inpt ncycle 10")
os.system("mod_map run.inpt ntstep 1")
#os.system("mod_map run.inpt rsdl_debug 1")

#os.system("mpiexec -np 3 tri_hp_petsc run.inpt -stop_for_debugger")
os.system("mpiexec -np 3 tri_hp_petsc run.inpt")


os.chdir("..")
os.system("opendiff Baseline/ Results/")
#os.system("open naca.pdf");