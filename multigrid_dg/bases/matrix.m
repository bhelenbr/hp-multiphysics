clear ms cv df leftl leftr rightl rightr
clear msc cvc dfc leftlc leftrc rightlc rightrc
clear bright bdright bleft bdleft dbasis lbasis rbasis

% CREATE MATRICES
% for i1=1:P+1
%     for i2=1:P+1
%         ms(i1,i2) = int(basis(i1)*basis(i2),-0.5,0.5);
%         cv(i1,i2) = -int(diff(basis(i1))*basis(i2),-0.5,0.5);
%         df(i1,i2) = int(diff(basis(i1))*diff(basis(i2)),-0.5,0.5);
%         leftl(i1,i2) =  subs(basis(i1),x,-0.5)*subs(basis(i2),x,0.5);
%         leftr(i1,i2) =  subs(basis(i1),x,-0.5)*subs(basis(i2),x,-0.5);
%         rightl(i1,i2) = subs(basis(i1),x,0.5)*subs(basis(i2),x,0.5);
%         rightr(i1,i2) = subs(basis(i1),x,0.5)*subs(basis(i2),x,-0.5);
%     end
%     bleft(i1,1) = subs(basis(i1),x,-0.5);
%     bright(i1,1) = subs(basis(i1),x,0.5);
%     bdleft(i1,1) = subs(diff(basis(i1)),x,-0.5);
%     bdright(i1,1) = subs(diff(basis(i1)),x,0.5);
% end

for i1=1:size(basis,2)
    dbasis(i1) = diff(basis(i1));
    lbasis(i1) = subs(basis(i1),x,-1);
    rbasis(i1) = subs(basis(i1),x,1);
end

% % CREATE MATRICES
% for i1=1:P+1
%     for i2 = 1:i1-1
%         ms(i1,i2) = int(basis(i1)*basis(i2),-1,1);
%         ms(i2,i1) = ms(i1,i2);
%         df(i1,i2) = int(dbasis(i1)*dbasis(i2),-1,1);
%         df(i2,i1) = df(i1,i2);
%         cv(i1,i2) = -int(dbasis(i1)*basis(i2),-1,1);
%         cv(i2,i1) = -cv(i1,i2) +rbasis(i1)*rbasis(i2)  -lbasis(i1)*lbasis(i2);
%     end
%     i2 = i1;
%     ms(i1,i2) = int(basis(i1)*basis(i2),-1,1);
%     df(i1,i2) = int(dbasis(i1)*dbasis(i2),-1,1);
%     cv(i1,i2) = -int(dbasis(i1)*basis(i2),-1,1);
% end

for i1=1:P+1
    for i2=1:P+1
        ms(i1,i2) = int(basis(i1)*basis(i2),-1,1);
        df(i1,i2) = int(dbasis(i1)*dbasis(i2),-1,1);
        cv(i1,i2) = -int(dbasis(i1)*basis(i2),-1,1);
        leftl(i1,i2) =  lbasis(i1)*rbasis(i2);
        leftr(i1,i2) =  lbasis(i1)*lbasis(i2);
        rightl(i1,i2) = rbasis(i1)*rbasis(i2);
        rightr(i1,i2) = rbasis(i1)*lbasis(i2);
    end
    bleft(i1,1) = lbasis(i1);
    bright(i1,1) = rbasis(i1);
    bdleft(i1,1) = subs(dbasis(i1),x,-1);
    bdright(i1,1) = subs(dbasis(i1),x,1);
end
ms = double(ms);
cv = double(cv);
df = double(df);
leftl = double(leftl);
leftr = double(leftr);
rightl = double(rightl);
rightr = double(rightr);

% CREATE COARSE MATRICES
msc = restrict*ms*restrict';
cvc = restrict*cv*restrict';
dfc = restrict*df*restrict';
leftlc = restrict*leftl*restrict';
leftrc = restrict*leftr*restrict';
rightlc = restrict*rightl*restrict';
rightrc = restrict*rightr*restrict';