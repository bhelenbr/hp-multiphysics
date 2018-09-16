#! /usr/bin/env python

# Runs a case with a facet point
# Flow is uniform flow
# All cases except symmetric perform similarly
# under_relaxation is needed to get log2p=0 to log2p=1 switch to work
# once it worked without it, but I never figured out what was different
# uncomment full_test to do the whole thing
# otherwise does one time step only
# For the short test if adapt is on, it will diverge at log2p = 2
# with adapt off it will converge (with no under_relaxation)
# symmetric has not made it all the way to the turning point
# Because v,p = 0, can get poor convergence because of numerical Jacobian
# Sensitive to eps_a in evaluating Jacobian when using dw = dw*eps_r +eps_a
# With new kinetic expression it stopped one time step before when using the old kinetic expression!

import sys
import os
import popen2
import numpy
import matplotlib.pyplot as plt
import subprocess
#import sr
import glob
import string
import math

os.chdir(os.path.dirname(sys.argv[0]))

FULL_TEST=1

HP="mpiexec -np 2 tri_hp_petsc run.inpt"

if not os.path.isdir("Results"):
	os.mkdir("Results")
os.chdir("Results")
os.system("rm *")

# copy input files into results directory
os.system("cp ../Inputs/* .")

# generate mesh and remove unnecessary data files
os.system("tri_mesh generate")
os.system("rm data*.grd")

# Various ways of running
#mod_map run.inpt b0_s3_one_sided 1
#mod_map run.inpt b0_s3_precondition 1

# This way never made it from log2p = 1 to log2p = 2
# mod_map run.inpt b0_s3_symmetric 1
# mod_map run.inpt b1_v1_hp_type hp_deformable_free_pnt
# mod_map run.inpt b1_v2_hp_type melt_facet_pt
	
# These are the steps to run
# Get a reasonable steady temperature field
# Before starting the melting 
os.system("cp run.inpt startup.inpt")
os.system("mod_map run.inpt b0_v1_hp_type plain")
os.system("mod_map run.inpt b0_v2_hp_type plain")
os.system("mod_map run.inpt b0_s3_hp_type inflow")
os.system("mod_map run.inpt b0_s3_type symbolic")
os.system("mod_map run.inpt b1_s3_hp_type dirichlet")
os.system("mod_map run.inpt b1_s3_type symbolic")
os.system("mod_map run.inpt b1_v1_hp_type plain")
os.system("mod_map run.inpt b1_v2_hp_type plain")
os.system("mod_map -c run.inpt mesh_movement")
os.system("mod_map run.inpt ntstep 1")
os.system('mod_map run.inpt dtinv 0.0')
os.system("mod_map run.inpt restart_interval 1")
# adapt has to be off so s3 doesn"t get messed up.
os.system("mod_map run.inpt adapt 0")
err = os.system(HP)
if err:
  print "exited from initial T solution"
  exit
  
RESTART=1

if FULL_TEST == 0:
	os.system("mv startup.inpt run.inpt")
	os.system("mod_map run.inpt adapt 0")
	os.system("mod_map run.inpt restart_interval 1")
	os.system("mod_map run.inpt ntstep 1")
	os.system("mod_map run.inpt restart " +str(RESTART))
	os.system('mod_map run.inpt dtinv_prev 0.0')
	# os.system("mod_map run.inpt debug_output 1")
	# os.system("mod_map run.inpt rsdl_debug 2")
	# os.system("mod_map run.inpt jac_debug 1")
	err = os.system(HP)
	if err:
  		print "exited from quick test with log2p = 0"
 		sys.exit()
	
	RESTART=RESTART+1
	
	print "log2p=1"
	os.system("mod_map run.inpt log2p 1")
	os.system("mod_map run.inpt restart " +str(RESTART))
	os.system("mod_map run.inpt -d dtinv_prev")
	err = os.system(HP)
	if err:
  		print "exited from quick test with log2p = 1"
 		sys.exit()
	
	RESTART=RESTART+1
	
	print "log2p=2"
	os.system("mod_map run.inpt log2p 2")
	os.system("mod_map run.inpt restart " +str(RESTART))
	err = os.system(HP)
	if err:
  		print "exited from quick test with log2p = 2"
 		sys.exit()
	
	
	os.chdir("..");
	os.system("opendiff Baseline/ Results/")
	sys.exit()

os.system("mv startup.inpt run.inpt")
# 100 steps with small time step
os.system("mod_map run.inpt ntstep 100")
os.system("mod_map run.inpt restart " +str(RESTART))
os.system("mod_map run.inpt restart_interval 100")
os.system('mod_map run.inpt dtinv_prev 0.0')
#os.system("mod_map run.inpt debug_output 1")
#os.system("mod_map run.inpt under_relaxation 0.66")
err = os.system(HP)
if err:
	print "exited from unsteady evolution 1"
	sys.exit()
 
RESTART=100
print "unsteady evolution 1 " +str(RESTART)


# More steps with larger time step
os.system("mod_map run.inpt restart " +str(RESTART))
os.system('mod_map run.inpt dtinv_prev "1e5*(dtinv1 +dtinv3)"')
os.system('mod_map run.inpt dtinv "1e4*(dtinv1+dtinv3)"')
err = os.system(HP)
if err:
	print "exited from unsteady evolution 2"
	sys.exit()

RESTART=RESTART+100
print "unsteady evolution 3 " +str(RESTART)
# 200

# More steps with larger time step
os.system("mod_map run.inpt restart " +str(RESTART))
os.system('mod_map run.inpt dtinv_prev "1e4*(dtinv1 +dtinv3)"')
os.system('mod_map run.inpt dtinv "1e3*(dtinv1+dtinv3)"');
err = os.system(HP)
if err:
	print "exited from unsteady evolution 3"
	sys.exit()

RESTART=RESTART+100
print "unsteady evolution 3 " +str(RESTART)
# 200


# More steps with larger time step
os.system("mod_map run.inpt restart " +str(RESTART))
os.system('mod_map run.inpt dtinv_prev "1e3*(dtinv1 +dtinv3)"')
os.system('mod_map run.inpt dtinv "1e2*(dtinv1+dtinv3)"')
err = os.system(HP)
if err:
	print "exited from unsteady evolution 4"
	sys.exit()

RESTART=RESTART+100
print "unsteady evolution 4 " +str(RESTART)
# 300

# More steps with larger time step
os.system("mod_map run.inpt restart " +str(RESTART))
os.system('mod_map run.inpt dtinv_prev "1e2*(dtinv1 +dtinv3)"')
os.system('mod_map run.inpt dtinv "1e1*(dtinv1+dtinv3)"')
err = os.system(HP)
if err:
	print "exited from unsteady evolution 5"
	sys.exit()

RESTART=RESTART+100
print "unsteady evolution 5 " +str(RESTART)
# 400

# Steady State
os.system("mod_map run.inpt restart " +str(RESTART))
os.system("mod_map run.inpt ntstep 1")
os.system("mod_map run.inpt restart_interval 1")
os.system("mod_map -u run.inpt under_relaxation")
os.system("mod_map -d run.inpt dtinv_prev")
os.system('mod_map run.inpt dtinv 0.0')
err = os.system(HP)
if err:
	print "exited from initial steady-state"
	sys.exit()

RESTART=RESTART+1
print "initial steady-state " +str(RESTART)
# 401

# MAKE K2DN_MAX infinite
os.system("mod_map run.inpt restart " +str(RESTART))
os.system("mod_map run.inpt b0_s3_K2Dn_max 1e10")
err = os.system(HP)
if err:
	print "exited from changing K2DN_max"
	sys.exit()

RESTART=RESTART+1
print "changing K2DN_max " +str(RESTART)
# 402

# Get high order solution (try to get 1step steady-state first)
os.system('mod_map run.inpt dtinv 0.0')
os.system("mod_map run.inpt log2p 1")
os.system("mod_map run.inpt restart " +str(RESTART))
os.system("mod_map run.inpt adapt 0")
os.system("mod_map run.inpt under_relaxation 0.25")
err = os.system(HP)
# if err:
#   print "exited from log2p = 1 steady state"
#   exit
# 
if err:
	os.system('mod_map run.inpt dtinv "1e4*(dtinv1 +dtinv3)"')
	os.system("mod_map run.inpt ntstep 20")
	os.system("mod_map -c run.inpt under_relaxation")
	err = os.system(HP)
	if err:
		print "exited from log2p = 1"
		sys.exit()
	
	RESTART=RESTART+20
	os.system('mod_map run.inpt dtinv 0.0')
	os.system("mod_map run.inpt restart " +str(RESTART))
	os.system("mod_map run.inpt ntstep 1")
	os.system("mod_map -u run.inpt under_relaxation")
	err = os.system(HP)
	if err:
		print "exited from log2p = 1"
		sys.exit()
	

RESTART=RESTART+1
print "log2p = 1  " +str(RESTART)

# Get high order solution (try to get 1step steady-state first)
os.system("mod_map run.inpt log2p 2")
os.system("mod_map run.inpt restart " +str(RESTART))
os.system("mod_map run.inpt adapt 0")
os.system('mod_map run.inpt dtinv 0.0')
os.system("mod_map run.inpt ntstep 1")
os.system("mod_map -u run.inpt under_relaxation")
err = os.system(HP)
# if err:
#   print "exited from log2p = 2 steady state"
#   exit
# 
if err:
	os.system('mod_map run.inpt dtinv "1e4*(dtinv1 +dtinv3)"')
	os.system("mod_map run.inpt ntstep 20")
	os.system("mod_map -c run.inpt under_relaxation")
	err = os.system(HP)
	if err:
		print "exited from log2p = 2"
		sys.exit()
	
	RESTART=RESTART+20
	os.system('mod_map run.inpt dtinv 0.0')
	os.system("mod_map run.inpt restart " +str(RESTART))
	os.system("mod_map run.inpt ntstep 1")
	os.system("mod_map -u run.inpt under_relaxation")
	err = os.system(HP)
	if err:
		print "exited from log2p = 1"
		sys.exit()
	

RESTART=RESTART+1
print "log2p = 2 " +str(RESTART)
os.system("mod_map run.inpt under_relaxation 1.0")

# Get steady-state
os.system('mod_map run.inpt dtinv 0.0')
os.system("mod_map run.inpt restart " +str(RESTART))
os.system("mod_map run.inpt ntstep 1")
err = os.system(HP)
if err:
	print "exited from log2p = 2 steady state"
	sys.exit()

RESTART=RESTART+1
print "log2p = 2 steady state " +str(RESTART)

# Turn on Mesh adaptation
os.system("mod_map -u run.inpt error_estimator")
os.system("mod_map run.inpt adapt 1")
os.system("mod_map run.inpt restart " +str(RESTART))
os.system("mod_map run.inpt ntstep 3")
err = os.system(HP)
if err:
	print "exited from initial mesh adaptation"
	sys.exit()

RESTART=RESTART+3
print "mesh adaptation " +str(RESTART)

################################################
# SX IS THE PULL SPEED
# NOW SLOWLY INCREASE SX
################################################
os.system("mod_map -c run.inpt under_relaxation")
os.system('mod_map run.inpt dtinv 0.0')
os.system("mod_map run.inpt ntstep 1")
os.system("delete_data.bash 2 0")

sx = numpy.zeros(200)
xle = numpy.zeros(200)

p0 = subprocess.Popen(["mod_map","-e","run.inpt","sx"], stdout=subprocess.PIPE)
p1 = subprocess.Popen(["cut","-d/","-f1"], stdin=p0.stdout, stdout=subprocess.PIPE)
(out, err) = p1.communicate()
sx[0] = float(out)
 
filename = "data" +str(RESTART) +"_b1.dat"
p0 = subprocess.Popen(["head","-2",filename], stdout=subprocess.PIPE)
p1 = subprocess.Popen(["tail","-1"], stdin=p0.stdout, stdout=subprocess.PIPE)
p2 = subprocess.Popen(["cut","-d"," ","-f","1"], stdin=p1.stdout, stdout=subprocess.PIPE)
(out, err) = p2.communicate()
xle[0] = float(out)

SXMAX=0.1
FACTOR=1.05
DSX=(FACTOR-1)*sx[0]
ATTEMPTS=0
MAXATTEMPTS=3

count = 1

while ( ATTEMPTS < MAXATTEMPTS and sx[count] < SXMAX):
	sx[count] = sx[count-1] +DSX
	print "SX " +str(sx[count])
	os.system('mod_map run.inpt sx "' +str(sx[count]) +'/(d0*tsi)"');
	os.system("mod_map run.inpt restart " +str(RESTART))
	err = os.system(HP)
	if err:
		os.system('rm core*')
		os.system('rm neg*')
		os.system('rm abort*')
		DSX /= 4.0;
		os.system('mod_map run.inpt extrapolate 0.25');
		ATTEMPTS+=1
	else:
		# TEST TO MAKE SURE CONVERGED
		p0 = subprocess.Popen(["grep","^[1-9]","out_b0.log"], stdout=subprocess.PIPE)
		p1 = subprocess.Popen(["tail","-1"], stdin=p0.stdout, stdout=subprocess.PIPE)
		p2 = subprocess.Popen(["cut","-d"," ","-f","3"], stdin=p1.stdout, stdout=subprocess.PIPE)
		(out, err) = p2.communicate()
		res = float(out)
		if res > 1.0e-5:
	  		DSX /= 4.0;
	  		os.system('mod_map run.inpt extrapolate 0.25');
	  		ATTEMPTS+=1
		else:
			RESTART+=1
			os.system('mod_map run.inpt extrapolate 1.0');
			filename = "data" +str(RESTART) +"_b1.dat"
			p0 = subprocess.Popen(["head","-2",filename], stdout=subprocess.PIPE)
			p1 = subprocess.Popen(["tail","-1"], stdin=p0.stdout, stdout=subprocess.PIPE)
			p2 = subprocess.Popen(["cut","-d"," ","-f","1"], stdin=p1.stdout, stdout=subprocess.PIPE)
			(out, err) = p2.communicate()
			xle[count] = float(out)
			if count > 2:
				alpha = numpy.zeros(3);
				alpha[2] = sx[count]/((xle[count]-xle[count-1])*(xle[count]-xle[count-2]))
				alpha[1] = sx[count-1]/((xle[count-1]-xle[count])*(xle[count-1]-xle[count-2]))
				alpha[0] = sx[count-2]/((xle[count-2]-xle[count])*(xle[count-2]-xle[count-1]))
				
				xleturn = alpha[0]*(-xle[count]-xle[count-1])+alpha[1]*(-xle[count-2]-xle[count])+alpha[2]*(-xle[count-1]-xle[count-2])
				xleturn /= -2*(alpha[0]+alpha[1]+alpha[2])
				sxturn = alpha[2]*(xleturn-xle[count-1])*(xleturn-xle[count-2]) +alpha[1]*(xleturn-xle[count])*(xleturn-xle[count-2]) +alpha[2]*(xleturn-xle[count-1])*(xleturn-xle[count-2])
				print str(sxturn)
				if (sxturn > sx[count] and FACTOR*sx[count] > sxturn):
					# DSX /= 4.0;
					# os.system('mod_map run.inpt extrapolate 0.25');
					print "would have changed factor"
			count += 1

plt.plot(sx[0:count],xle[0:count],'r-x')
plt.xlabel('pull speed / [m/s]')
plt.ylabel('leading edge position')
plt.savefig("turning.pdf")

data = numpy.zeros([200,2])
data[:,1] = sx
data[:,2] = xle
numpy.savetxt("turning.dat",data)


os.system('rm core* abort* net* rstrt*.nc')
