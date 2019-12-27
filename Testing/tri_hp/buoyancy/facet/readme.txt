How to interpret files:

dataXXX_b0.dat:  Tecplot file for timestep XXX.  This is the tecplot file for the liquid.  There are 6 variables,  x, y, u, v, T, and p.   

dataXXX_b1.dat:  Tecplot file for timestep XXX.  This is the tecplot file for the solid.  There are 3 variables, x, y, and T.

dataXXX_b0_s3.dat:  1D data give solution values and gradients on interface.   Variables are "S", "X", "Y", "DXDT", "DYDT", "V0", "DVTANG0", "DVNORM0", "V1", "DVTANG1", "DVNORM1", "V2", "DVTANG2", "DVNORM2", "V3", "DVTANG3", "DVNORM3".   S is the arclength measured from the exit point (not so useful for this problem) X, Y are position, DXDT & DYDT are velocities of the interface at that point.   The rest are the values and tangential and normal derivatives of u, v, T, and P.

dataXXX_b0_s4.dat:  Same format as above except this is the surface of the liquid.

dataXXX_b1_s3.dat:  Same format as above except this is for the interface on the solid side.  There is only one variable and its derivatives for temperature.

dataXXX_b1_s7.dat:  Same format as above except this is for the top surface of the solid.

dataXXXkinetics_b0_s3: "S", "X", "Y", "K", "Delta T", and "sin(theta)"

Temperatures are nondimensional by T_melt.   Distance by depth of the crucible.  To get heat fluxes you have to multiply Gradients by the thermal conductivity and divide by the crucible depth.