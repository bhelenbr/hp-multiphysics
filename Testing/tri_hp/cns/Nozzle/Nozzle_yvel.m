%Initial condition for v

syms M gam RTo Po
rhoo = Po/RTo;

Ar = (1/M)*((2+(gam-1)*M^2)/(gam+1))^((gam+1)/(2*(gam-1)));
dArdM = diff(Ar,M);

RT = RTo/(1+0.5*(gam-1)*M^2);
P = Po/((1+0.5*(gam-1)*M^2)^(gam/(gam-1)));

c = sqrt(gam*RT);
u = M*c;
dudM = diff(u,M);

dPdM = diff(P,M);
dRTdM = diff(RT,M);
dcdM = diff(c,M);
