syms x real;

% FOR DUBINER BASIS
basis(1) = (1-x)/2;
basis(2) = (1+x)/2;
for n = 3:P+1
    basis(n) = simplify((1-x)/2*(1+x)/2*maple('JacobiP',n-3, 1, 1, x));
end

% RESTRICTION FOR HEIRARCHICAL SYSTEM
restrict = eye(PC+1,P+1);

% ROTATE LAST MODE
btemp = basis(2);
basis(2:P) = basis(3:P+1);
basis(P+1) = btemp;
restrict = mrotate(restrict);

matrix;