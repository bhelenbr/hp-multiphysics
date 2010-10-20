clear basis2d ms2d cvx cvy dfx dfy dyx dxy ll2d lr2d rl2d rr2d;
clear bb2d bt2d tb2d tt2d brdr bldl bldr brdl btdt bbdb bbdt btdb;
clear restrict2d;
        
% ORDERING OF POLYNOMIALS IS IN STRIPS IN X
nbasis = (P+1);
nbasisc = (PC+1);
basis2d = basis;

ms2d = ms;
cvx = cv;
cvy = zeros(size(cvx));
dfx = df;
dfy = zeros(size(dfx));
dxy = zeros(size(cvx));
dyx = zeros(size(cvx));

ll2d = leftl;
lr2d = leftr;
rl2d = rightl;
rr2d = rightr;

bb2d = zeros(size(ll2d));
bb2d = zeros(size(ll2d));
bt2d = zeros(size(ll2d));
tb2d = zeros(size(ll2d));
tt2d = zeros(size(ll2d));

brdr = bright*bdright';
bldl = bleft*bdleft';
bldr = bleft*bdright';
brdl = bright*bdleft';      

btdt = zeros(size(ll2d));
bbdb = zeros(size(ll2d));
bbdt = zeros(size(ll2d));
btdb = zeros(size(ll2d));

restrict2d= restrict;


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