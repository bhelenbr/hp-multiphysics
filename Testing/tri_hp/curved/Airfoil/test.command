#! /usr/bin/env python

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
p0 = subprocess.Popen("echo ${PWD%/*/*}/bin/:", stdout=subprocess.PIPE,shell=True)
(BINDIR, err) = p0.communicate()
os.environ['PATH'] = BINDIR[:-1].decode('ascii') + os.environ['PATH']

if not os.path.isdir("Results"):
	os.mkdir("Results")
os.chdir("Results")
os.system("rm *")

# copy input files into results directory
os.system("cp ../Inputs/* .")

nlayers = 10
# generate mesh and remove unnecessary data files
offset = 0.0
#os.system("spline -m 0.25,0.0 -r -10.0 -s 3.0 -o" +str(offset)+ " -i spoints.dat naca.spl > interp.dat");
os.system("spline -m 0.0,0.0 -r -0.0 -s 1.0 -o " +str(offset)+ " naca.spl 60 3.0 >> interp.dat");

dy = -1e-3
offset = dy
for x in range(nlayers):
	#os.system("spline -m 0.25,0.0 -r -10.0 -s 3.0 -o" +str(offset)+ " -i spoints.dat naca.spl > interp.dat");
	os.system("spline -m 0.0,0.0 -r -0.0 -s 1.0 -o " +str(offset)+ " naca.spl 60 3.0 >> interp.dat");
	dy *= 1.3
	offset += dy


import sys
import numpy as np
import matplotlib.pyplot as plt
import math

s, x, y, tx, ty, curvx, curvy = np.loadtxt("interp.dat", delimiter=' ', unpack=True)


plt.plot(x,y,'r-x')
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')
plt.savefig("naca.pdf")

# output .d files
npoints = int(s.size/(nlayers+1))
f = open('blayer.d','w')
f.write(str(s.size-1+2*nlayers+1)+'\n')

count = 0
# airfoil surface
#skip trailing edge point
for i in range(npoints-1):
	f.write('{0:3d}: {1:.10f} {2:.10f} 0.01 0\n'.format(count,x[i],y[i]))
	count += 1

# Outer surface 
for i in range(s.size-npoints,s.size):
	f.write('{0:3d}: {1:.10f} {2:.10f} {3:.3f} 0\n'.format(count,x[i],y[i],math.sqrt(math.pow(x[i]-x[i-1],2)+math.pow(y[i]-y[i-1],2))))
	count += 1
	
#add one column of wake points (top)
for i in range(0,nlayers):
	dx = y[(i+1)*npoints+npoints-1] -y[i*npoints+npoints-1]
	f.write('{0:3d}: {1:.10f} {2:.10f} {3:.3f} 0\n'.format(count,x[i*npoints+npoints-1]+dx,y[i*npoints+npoints-1],dx))
	count += 1
i = nlayers
f.write('{0:3d}: {1:.10f} {2:.10f} {3:.3f} 1\n'.format(count,x[i*npoints+npoints-1]+dx,y[i*npoints+npoints-1],dx))
count += 1
	
#add one column of wake points (bottom)
for i in range(1,nlayers):
	dx = y[i*npoints] -y[(i+1)*npoints]
	f.write('{0:3d}: {1:.10f} {2:.10f} {3:.3f} 0\n'.format(count,x[i*npoints]+dx,y[i*npoints],dx))
	count += 1
i = nlayers
f.write('{0:3d}: {1:.10f} {2:.10f} {3:.3f} 2\n'.format(count,x[i*npoints]+dx,y[i*npoints],dx))
count += 1

#write interior points
for i in range(npoints,s.size-npoints):
	f.write('{0:3d}: {1:.10f} {2:.10f} 0.01 0\n'.format(count,x[i],y[i]))
	count += 1
	
f.write('{0:d}\n'.format(2*(npoints-1)+2*nlayers+2))
# Airfoil surface
count = 0
for i in range(npoints-2):
	f.write('{0:d}: {1:d} {2:d} 1\n'.format(count,i,i+1))
	count += 1
f.write('{0:d}: {1:d} {2:d} 1\n'.format(count,npoints-2,0))
count += 1

# Outer boundary
for i in range(npoints-1):
	f.write('{0:d}: {1:d} {2:d} 2\n'.format(count,2*npoints-2-i,2*npoints-3-i))
	count += 1
# connection to wake point top
f.write('{0:d}: {1:d} {2:d} 2\n'.format(count,2*npoints-1+nlayers,2*npoints-2))
count += 1
# connection to wake point bot
f.write('{0:d}: {1:d} {2:d} 2\n'.format(count,2*npoints-1-npoints,2*npoints-1+2*nlayers))
count += 1

# Top trailing edge
for i in range(nlayers):
	f.write('{0:d}: {1:d} {2:d} 3\n'.format(count,2*npoints-1+i,2*npoints-1+i+1))
	count += 1

# Bottom trailing edge
f.write('{0:d}: {1:d} {2:d} 3\n'.format(count,2*npoints-1+nlayers+1,2*npoints-1))
count += 1
for i in range(1,nlayers):
	f.write('{0:d}: {1:d} {2:d} 3\n'.format(count,2*npoints-1+nlayers+i+1,2*npoints-1+nlayers+i))
	count += 1

f.close()


f = open('domain.d','w')
f.write(str(npoints +2*nlayers+1 +4)+'\n')

count = 0
# Outer surface 
for i in range(s.size-npoints,s.size):
	f.write('{0:3d}: {1:.10f} {2:.10f} {3:.3f} 0\n'.format(count,x[i],y[i],math.sqrt(math.pow(x[i]-x[i-1],2)+math.pow(y[i]-y[i-1],2))))
	count += 1
	
#add one column of wake points (top)
for i in range(nlayers):
	dx = y[(i+1)*npoints+npoints-1] -y[i*npoints+npoints-1]
	f.write('{0:3d}: {1:.10f} {2:.10f} {3:.3f} 0\n'.format(count,x[i*npoints+npoints-1]+dx,y[i*npoints+npoints-1],dx))
	count += 1
i = nlayers
f.write('{0:3d}: {1:.10f} {2:.10f} {3:.3f} 1\n'.format(count,x[i*npoints+npoints-1]+dx,y[i*npoints+npoints-1],dx))
count += 1
	
#add one column of wake points (bottom)
for i in range(1,nlayers):
	dx = y[i*npoints] -y[(i+1)*npoints]
	f.write('{0:3d}: {1:.10f} {2:.10f} {3:.3f} 0\n'.format(count,x[i*npoints]+dx,y[i*npoints],dx))
	count += 1
i = nlayers
f.write('{0:3d}: {1:.10f} {2:.10f} {3:.3f} 2\n'.format(count,x[i*npoints]+dx,y[i*npoints],dx))
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
	
f.write('{0:d}\n'.format(npoints+1+2*nlayers+4))
# Outer boundary layer surface
count = 0
for i in range(npoints-1):
	f.write('{0:d}: {1:d} {2:d} 2\n'.format(count,i,i+1))
	count += 1
# connection to wake point top
f.write('{0:d}: {1:d} {2:d} 2\n'.format(count,npoints-1,npoints-1+nlayers+1))
count += 1
# connection to wake point bot
f.write('{0:d}: {1:d} {2:d} 2\n'.format(count,npoints+2*nlayers,0))
count += 1

# Top trailing edge
for i in range(nlayers):
	f.write('{0:d}: {1:d} {2:d} 3\n'.format(count,npoints+i+1,npoints+i))
	count += 1

# Bottom trailing edge
f.write('{0:d}: {1:d} {2:d} 3\n'.format(count,npoints,npoints+nlayers+1))
count += 1
for i in range(1,nlayers):
	f.write('{0:d}: {1:d} {2:d} 3\n'.format(count,npoints+nlayers+i,npoints+nlayers+i+1))
	count += 1
	
# Outer boundary
f.write('{0:d}: {1:d} {2:d} 5\n'.format(count,npoints+2*nlayers+1,npoints+2*nlayers+2))
count += 1
f.write('{0:d}: {1:d} {2:d} 6\n'.format(count,npoints+2*nlayers+2,npoints+2*nlayers+3))
count += 1
f.write('{0:d}: {1:d} {2:d} 7\n'.format(count,npoints+2*nlayers+3,npoints+2*nlayers+4))
count += 1
f.write('{0:d}: {1:d} {2:d} 8\n'.format(count,npoints+2*nlayers+4,npoints+2*nlayers+1))
count += 1

f.close()

os.system("tri_mesh generate.inpt")
os.system("mpiexec -np 1 tri_hp_petsc deform.inpt")
os.system("mv rstrt1_d0_b0.nc curvatures_b0.nc")
os.system("rm data*.dat")
os.system("mpiexec -np 2 tri_hp_petsc run.inpt")

#os.system("mpiexec -np 2 tri_hp_petsc run.inpt -stop_for_debugger")
