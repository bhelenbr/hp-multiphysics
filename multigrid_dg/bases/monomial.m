 syms x real;
 
 % FOR MONOMIAL BASIS
 for n =1:P+1
     basis(n) = x^(n-1);
 end
 
 % RESTRICTION FOR HEIRARCHICAL SYSTEM
 restrict = eye(PC+1,P+1);
 
 matrix;

% WITHOUT USING SYMBOLIC TOOLBOX
clear ms cv df leftl leftr rightl rightr
clear msc cvc dfc leftlc leftrc rightlc rightrc
clear bright bdright bleft bdleft dbasis lbasis rbasis

for i1=1:P+1
    lbasis(i1) = (-1)^(i1-1);
    rbasis(i1) = (1)^(i1-1);
end

% CREATE MATRICES
for i1=1:P+1
    for i2=1:P+1
        ms(i1,i2) = 1/(i1+i2-1)*(1 - (-1)^(i1+i2-1));
        if (i1 == 1 || i2 == 1)
            df(i1,i2) = 0;
        else
            df(i1,i2) = (i1-1)*(i2-1)/(i1+i2-3)*(1 -(-1)^(i1+i2-3));
        end
        if (i1 == 1)
            cv(i1,i2) = 0;
        else
            cv(i1,i2) = -(i1-1)/(i1+i2-2)*(1 -(-1)^(i1+i2-2));
        end
        leftl(i1,i2) =  lbasis(i1)*rbasis(i2);
        leftr(i1,i2) =  lbasis(i1)*lbasis(i2);
        rightl(i1,i2) = rbasis(i1)*rbasis(i2);
        rightr(i1,i2) = rbasis(i1)*lbasis(i2);
    end
    bleft(i1,1) = lbasis(i1);
    bright(i1,1) = rbasis(i1);
    bdleft(i1,1) = (i1-1)*(-1)^(i1-2);
    bdright(i1,1) = (i1-1)*(1)^(i1-2);
end

restrict = eye(PC+1,P+1);

% CREATE COARSE MATRICES
msc = restrict*ms*restrict';
cvc = restrict*cv*restrict';
dfc = restrict*df*restrict';
leftlc = restrict*leftl*restrict';
leftrc = restrict*leftr*restrict';
rightlc = restrict*rightl*restrict';
rightrc = restrict*rightr*restrict';

 
