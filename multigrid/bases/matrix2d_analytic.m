clear ms2d cvx cvy dfx dfy dyx dxy ll2d lr2d rl2d rr2d;
clear bb2d bt2d tb2d tt2d brdr bldl bldr brdl btdt bbdb bbdt btdb;

nbasis = size(basis2d,2);

for i1=1:nbasis
    for i2=1:nbasis
        % CHECK MATRICES
        ms2d(i1,i2) = int(int(basis2d(i1)*basis2d(i2),x,-1,1),y,-1,1);
        cvx(i1,i2) = -int(int(diff(basis2d(i1),x)*basis2d(i2),x,-1,1),y,-1,1);
        cvy(i1,i2) = -int(int(diff(basis2d(i1),y)*basis2d(i2),x,-1,1),y,-1,1);

        dfx(i1,i2) = int(int(diff(basis2d(i1),x)*diff(basis2d(i2),x),x,-1,1),y,-1,1);
        dfy(i1,i2) = int(int(diff(basis2d(i1),y)*diff(basis2d(i2),y),x,-1,1),y,-1,1);
        
        ll2d(i1,i2) = int(subs(basis2d(i1),x,-1)*subs(basis2d(i2),x,1),y,-1,1);
        lr2d(i1,i2) = int(subs(basis2d(i1),x,-1)*subs(basis2d(i2),x,-1),y,-1,1);
        rl2d(i1,i2) = int(subs(basis2d(i1),x,1)*subs(basis2d(i2),x,1),y,-1,1);
        rr2d(i1,i2) = int(subs(basis2d(i1),x,1)*subs(basis2d(i2),x,-1),y,-1,1);  
        
        bb2d(i1,i2) = int(subs(basis2d(i1),y,-1)*subs(basis2d(i2),y,1),x,-1,1);
        bt2d(i1,i2) = int(subs(basis2d(i1),y,-1)*subs(basis2d(i2),y,-1),x,-1,1);
        tb2d(i1,i2) = int(subs(basis2d(i1),y,1)*subs(basis2d(i2),y,1),x,-1,1);
        tt2d(i1,i2) = int(subs(basis2d(i1),y,1)*subs(basis2d(i2),y,-1),x,-1,1);
        
        brdr(i1,i2) = int(subs(basis2d(i1),x,1)*subs(diff(basis2d(i2),x),x,1),y,-1,1);
        bldl(i1,i2) = int(subs(basis2d(i1),x,-1)*subs(diff(basis2d(i2),x),x,-1),y,-1,1);
        bldr(i1,i2) = int(subs(basis2d(i1),x,-1)*subs(diff(basis2d(i2),x),x,1),y,-1,1);
        brdl(i1,i2) = int(subs(basis2d(i1),x,1)*subs(diff(basis2d(i2),x),x,-1),y,-1,1);
        
        btdt(i1,i2) = int(subs(basis2d(i1),y,1)*subs(diff(basis2d(i2),y),y,1),x,-1,1);
        bbdb(i1,i2) = int(subs(basis2d(i1),y,-1)*subs(diff(basis2d(i2),y),y,-1),x,-1,1);
        bbdt(i1,i2) = int(subs(basis2d(i1),y,-1)*subs(diff(basis2d(i2),y),y,1),x,-1,1);
        btdb(i1,i2) = int(subs(basis2d(i1),y,1)*subs(diff(basis2d(i2),y),y,-1),x,-1,1);
    end
end

ms2d = double(ms2d);
cvx = double(cvx);
cvy = double(cvy);

dfx = double(dfx);
dfy = double(dfy);

ll2d = double(ll2d);
lr2d = double(lr2d);
rl2d = double(rl2d);
rr2d = double(rr2d);

bb2d = double(bb2d);
bt2d = double(bt2d);
tb2d = double(tb2d);
tt2d = double(tt2d);

brdr = double(brdr);
bldl = double(bldl);
bldr = double(bldr);
brdl = double(brdl);

btdt = double(btdt);
bbdb = double(bbdb);
bbdt = double(bbdt);
btdb = double(btdb);

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