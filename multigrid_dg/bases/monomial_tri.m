clear basis2d restrict2d
syms x y real

indx = 1;
for m=0:P
    for n=0:m
        basis2d(indx) = x^(m-n)*y^(n);
        indx = indx +1;
    end
end

restrict2d = eye((PC+1)*(PC+2)/2,(P+1)*(P+2)/2);
    