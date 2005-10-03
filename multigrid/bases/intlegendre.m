syms x s real;

clear tempb
% FOR INTEGRATED LEGENDRE BASIS
basis(1) = sym(1);
for n = 1:P
    tempb = simplify(maple('LegendreP',n-1, s));
    basis(n+1) = int(tempb,0,x);
end

% RESTRICTION FOR HEIRARCHICAL SYSTEM
restrict = eye(PC+1,P+1);

matrix;