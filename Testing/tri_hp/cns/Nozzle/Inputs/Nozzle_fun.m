function F = Nozzle_fun(x)
A = x(1);
B = x(2);
C = x(3);
D = x(4);

F(1) = 2.5-A-C;
F(2) = -1+A-C*exp(-5*B)*cos(5*D-pi);
F(3) = -1.5+A-C*exp(-10*B)*cos(10*D-pi);
F(4) = -B*C*exp(-5*B)*cos(5*D)-C*D*exp(-5*B)*sin(5*D);