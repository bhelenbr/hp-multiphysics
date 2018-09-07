lambda = 2;

eta = 0:0.01:1

theta = 1 - erf(eta*lambda)/erf(lambda);

plot(eta,theta)
xlabel('\eta');
ylabel('\Theta');

myexport('analytic');

