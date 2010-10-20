clear ms cv df leftl leftr rightl rightr
clear msc cvc dfc leftlc leftrc rightlc rightrc

% CREATE MASS MATRIX
for i1=1:TM
    for i2=1:TM
        ms(i1,i2) = int(int(basis(i1)*basis(i2),s,-1,-r),r,-1,1);
    end
end
% ms = double(ms);
% cv = double(cv);
% df = double(df);
% 
% % CREATE COARSE MATRICES
% msc = restrict*ms*restrict';
% cvc = restrict*cv*restrict';
% dfc = restrict*df*restrict';
