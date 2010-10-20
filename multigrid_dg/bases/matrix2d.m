clear basis2d ms2d cvx cvy dfx dfy dyx dxy ll2d lr2d rl2d rr2d;
clear bb2d bt2d tb2d tt2d brdr bldl bldr brdl btdt bbdb bbdt btdb;
clear restrict2d;
        
% ORDERING OF POLYNOMIALS IS IN STRIPS IN X
nbasis = (P+1)*(P+1);
nbasisc = (PC+1)*(PC+1);

syms y

for i1=1:nbasis
    indx1 = rem(i1-1,P+1) +1;
    indy1 = floor((i1-1)/(P+1)) +1;
    basis2d(i1) = basis(indx1)*subs(basis(indy1),x,y);

	for i2=1:nbasis
        indx2 = rem(i2-1,P+1) +1;
        indy2 = floor((i2-1)/(P+1)) +1;
        ms2d(i1,i2) = ms(indx1,indx2)*ms(indy1,indy2);
        cvx(i1,i2) = cv(indx1,indx2)*ms(indy1,indy2);
        cvy(i1,i2) = ms(indx1,indx2)*cv(indy1,indy2);
        dfx(i1,i2) = df(indx1,indx2)*ms(indy1,indy2);
        dfy(i1,i2) = ms(indx1,indx2)*df(indy1,indy2);
        dxy(i1,i2) = cv(indx1,indx2)*cv(indy2,indy1);
        dyx(i1,i2) = cv(indx2,indx1)*cv(indy1,indy2);
        
        ll2d(i1,i2) = leftl(indx1,indx2)*ms(indy1,indy2);
        lr2d(i1,i2) = leftr(indx1,indx2)*ms(indy1,indy2);
        rl2d(i1,i2) = rightl(indx1,indx2)*ms(indy1,indy2);
        rr2d(i1,i2) = rightr(indx1,indx2)*ms(indy1,indy2);  
        
        bb2d(i1,i2) = leftl(indy1,indy2)*ms(indx1,indx2);
        bt2d(i1,i2) = leftr(indy1,indy2)*ms(indx1,indx2);
        tb2d(i1,i2) = rightl(indy1,indy2)*ms(indx1,indx2);
        tt2d(i1,i2) = rightr(indy1,indy2)*ms(indx1,indx2);
        
        brdr(i1,i2) = bright(indx1)*bdright(indx2)*ms(indy1,indy2);
        bldl(i1,i2) = bleft(indx1)*bdleft(indx2)*ms(indy1,indy2);
        bldr(i1,i2) = bleft(indx1)*bdright(indx2)*ms(indy1,indy2);
        brdl(i1,i2) = bright(indx1)*bdleft(indx2)*ms(indy1,indy2);
        
        btdt(i1,i2) = bright(indy1)*bdright(indy2)*ms(indx1,indx2);
        bbdb(i1,i2) = bleft(indy1)*bdleft(indy2)*ms(indx1,indx2);
        bbdt(i1,i2) = bleft(indy1)*bdright(indy2)*ms(indx1,indx2);
        btdb(i1,i2) = bright(indy1)*bdleft(indy2)*ms(indx1,indx2);
    end

    for i2 = 1:nbasisc;
        indx2 = rem(i2-1,PC+1) +1;
        indy2 = floor((i2-1)/(PC+1)) +1;
        restrict2d(i2,i1) = restrict(indx2,indx1)*restrict(indy2,indy1);
    end
end

% CREATE COARSE MATRICES
ms2dc = restrict2d*ms2d*restrict2d';
cvxc = restrict2d*cvx*restrict2d';
cvyc = restrict2d*cvy*restrict2d';
dfxc = restrict2d*dfx*restrict2d';
dfyc = restrict2d*dfy*restrict2d';

ll2dc = restrict2d*ll2d*restrict2d';
lr2dc = restrict2d*lr2d*restrict2d';
rl2dc = restrict2d*rl2d*restrict2d';
rr2dc = restrict2d*rr2d*restrict2d';  

bb2dc = restrict2d*bb2d*restrict2d';
bt2dc = restrict2d*bt2d*restrict2d';
tb2dc = restrict2d*tb2d*restrict2d';
tt2dc = restrict2d*tt2d*restrict2d';  

% % THIS IS TO TEST MATRIX2D 
% clear basis2d restrict2d;
% syms y;
% nbasis = (P+1)*(P+1);
% nbasisc = (PC+1)*(PC+1);
% for i1=1:nbasis
%     indx1 = rem(i1-1,P+1) +1;
%     indy1 = floor((i1-1)/(P+1)) +1;
%     basis2d(i1) = basis(indx1)*subs(basis(indy1),x,y);
%     for i2 = 1:nbasisc;
%         indx2 = rem(i2-1,PC+1) +1;
%         indy2 = floor((i2-1)/(PC+1)) +1;
%         restrict2d(i2,i1) = restrict(indx2,indx1)*restrict(indy2,indy1);
%     end
% end
% matrix2d_analytic;

