% SYMMETRIC GAUSS-SEIDEL ANALYSIS
syms k kl kll kr krr
syms m ml mll mr mrr

NEL = 8

for j=1:NEL
    K(j,j) = k;
    L(j,j) = m;
    U(j,j) = m;
end

for j=1:NEL-1
    K(j+1,j) = kl;
    K(j,j+1) = kr;
    L(j+1,j) = kl;
    U(j,j+1) = kr;
end

K(1,NEL) = kl;
K(NEL,1) = kr;
L(1,NEL) = kl;
U(NEL,1) = kr;

Rinv = (inv(U)+inv(L) -inv(U)*K*inv(L));
R = inv(Rinv)