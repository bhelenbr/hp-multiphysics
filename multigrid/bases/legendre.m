syms x real;

% FOR LEGENDRE BASIS
for n = 1:P+1
    basis(n) = simplify(maple('LegendreP',n-1, x));
end

% RESTRICTION FOR HEIRARCHICAL SYSTEM
restrict = eye(PC+1,P+1);

matrix;







