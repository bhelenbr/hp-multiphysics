fun = @Nozzle_fun;
z0 = [11/8 log(3)/5 9/8 0.2];
z = fsolve(fun,z0);

i = 1;
lam = 0.1*log(2.5/1.5);
for x = 0:0.1:10
    area(i) = z(1)-z(3)*exp(-z(2)*x)*cos(z(4)*x-pi);
    i = i+1;
end

plot(0:0.1:10,area)
hold on
plot(0:0.1:10,-area)

%Height at shock
syms x0
A1 = z(1)-z(3)*exp(-z(2)*x0)*cos(z(4)*x0-pi);
y = double(subs(A1,x0,7))

%Derivative
syms A B C D
A1 = A-C*exp(-B*x0)*cos(D*x0-pi);
diff(A1,x0);