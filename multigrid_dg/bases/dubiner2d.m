clear all;
P = input('P\n');
syms x n r s real;

% FOR DUBINER TRIANGULAR BASIS
basisx(1) = sym('1');
basisx(2) = (1-x)/2;
basisx(3) = (1+x)/2;
for m = 3:P+1
    basisx(m+1) = simplify((1-x)/2*(1+x)/2*maple('JacobiP',m-3, 2, 2, x));
end
basisx = simplify(basisx);

for i1=2:P+2
    for i2=2:P+2
        ms1d(i1-1,i2-1) = int(basisx(i1)*basisx(i2),x,-1,1);
    end
end

basisn(1) = (1+n)/2;
basisn(2) = (1-n)/2;
basisn(3) = (1-n)/2;

% SIDE 0 
for m = 2:P
  basisn(m+2) = ((1-n)/2)^m;
end

% SIDE 1  
for m=0:P-2
  basisn(P+3+m) = (1-n)/2*(1+n)/2*maple('JacobiP',m, 2, 2, n);
end

% SIDE 2  
for m=0:P-2
	basisn(2*P +2+m) = -1^m*(1-n)/2*(1+n)/2*maple('JacobiP',m, 2, 2, n);
end

indx = 2*P+2+P-2;
% INTERIOR MODES
for m = 2:P-1     
	for k = 0: P-1-m
        indx = indx+1;
		basisn(indx) = ((1-n)/2)^m*(1+n)/2*maple('JacobiP',k, 2*m+1, 2, n);
	end
end
basisn = simplify(basisn);

basisxn(1) = basisx(1)*basisn(1);
basisxn(2) = basisx(2)*basisn(2);
basisxn(3) = basisx(3)*basisn(3);

indx = 3;
% SIDE 1
S1 = indx+1;
for m=0:P-2
    indx = indx+1;
    basisxn(indx) = basisx(m+4)*basisn(m+4);
end

% SIDE 2
S2 = indx+1;
for m=0:P-2
    indx = indx+1;
    basisxn(indx) = basisx(3)*basisn(P+3+m);
end

% SIDE 3
S3 = indx+1;
for m=0:P-2
    indx = indx+1;
    basisxn(indx) = basisx(2)*basisn(2*P +2+m);
end

% INTERIOR MODES
BM = indx;
for m = 2:P-1     
	for k = 0: P-1-m
        indx = indx+1;
		basisxn(indx) = basisx(m+2)*basisn(indx);
	end
end
basisxn = simplify(basisxn);
basis = simplify(subs(basisxn,{x,n},{(2*(1+r)/(1-s) -1),s}));

for m=1:41
    xplotm(m,:) = -1:1/20:1;
    nplotm(m,:) = -1:1/20:1;
end
nplotm = nplotm';
rplotm = (xplotm + ones(size(xplotm))).*(ones(size(nplotm)) -nplotm)/2  -ones(size(nplotm));
splotm = nplotm;

% for m=1:indx
% 	figure;
% 	axis equal;
% 	zplotm = double(subs(basis(m),{r,s},{rplotm,splotm}));
% 	mesh(rplotm,splotm,zplotm)
% 	ezsurf((x+1)*(1-n)/2 -1,n,basisxn(m),[-1,1]);
% end

TM = indx;
matrix_tri;

% RENORMALIZE INTERIOR MODES 
for i=BM+1:TM
    basis(i) = basis(i)/sqrt(ms(i,i));
end

matrix_tri;

% lumped = basis;
% 
% % FIND MASS LUMPNG VERTEX MODE
% % CONSTRAINTS
% clm = [2,3];
% clm = [clm,S1:S1+P-3,S2:S2+P-3,S3:S3+P-3];
% indx = BM;
% for m = 2:P-1     
% 	for k = 0: P-2-m
%         indx = indx+1;
%         clm = [clm,indx];
% 	end
%     indx = indx +1;
% end
% 
% % FREEDOM
% row = [S2:S2+P-2,S3:S3+P-2];
% row = [row,BM+1:TM];
% 
% a = ms(row,clm);
% a = a';
% b = ms(1,clm)';
% 
% temp = a\b;
% ifmv(1,row) = temp;
% 
% lumped(1) = basis(1) -ifmv(1,:)*(basis');
% 
% % FIND SIDE 0 MODES!!!!
% 
% for sm=0:P-3
% 	% CONSTRAINTS
% 	clm = [S2+sm:S2+P-3,S3+sm:S3+P-3];
% 	skip = BM;
%     for m = 2:sm+1
%         for k = 0: P-1-m
%             skip = skip+1;
% 		end
% 	end
%     indx = skip;
% 	for m = sm+2:P-1     
% 		for k = 0: P-2-m
%             indx = indx+1;
%             clm = [clm,indx];
% 		end
%         indx = indx +1;
% 	end
% 	
% 	row = [S1+sm+1:S1+P-2];
% 	row = [row,skip+1:TM];
% 	
% 	a = ms(row,clm);
% 	a = a';
% 	b = ms(S1+sm,clm)';
% 	
% 	temp = a\b;
% 	ifmv(S1+sm,row) = temp;
% 	
% 	lumped(S1+sm) = basis(S1+sm) -ifmv(S1+sm,:)*(basis');
% end
% 
% figure;
% axis equal;
% zplotm = double(subs(lumped(1),{r,s},{rplotm,splotm}));
% contourf(rplotm,splotm,zplotm);
% box off;
% 
% v1d = subs(lumped(1),r,-1);
% %maple('evalf',solve(v1d))
% figure;
% ezplot(v1d,[-1,1]);
% 
% 
% 
% % TO ROTATE MODES
% %lumped(2) = subs(lumped(1),{r,s},{s,r+s});
% syms u v
% lumped(2) = subs(lumped(1),{r,s},{u,v});
% lumped(2) = subs(lumped(2),{u,v},{s,-1-r-s});
% lumped(3) = subs(lumped(2),{r,s},{u,v});
% lumped(3) = subs(lumped(3),{u,v},{s,-1-r-s});
% 
% for sm=0:P-2
%     lumped(S2+sm) = subs(lumped(S1+sm),{r,s},{u,v});
%     lumped(S2+sm) = subs(lumped(S2+sm),{u,v},{s,-1-r-s});
%     lumped(S3+sm) = subs(lumped(S2+sm),{r,s},{u,v});
%     lumped(S3+sm) = subs(lumped(S3+sm),{u,v},{s,-1-r-s});
% end
% 
% figure;
% axis equal;
% zplotm = double(subs(lumped(3),{r,s},{rplotm,splotm}));
% contourf(rplotm,splotm,zplotm);
% box off;
% 
% for sm=0:P-2
% 	figure;
% 	axis equal;
% 	zplotm = double(subs(lumped(S3+sm),{r,s},{rplotm,splotm}));
%     box off;
% 	contourf(rplotm,splotm,zplotm);
% end 

% FIND ORTHOGONAL TRIANGULARLY SYMMETRIC INTERIOR MODES
% FIRST ROTATE BASIS
syms u v
int_rot(1:BM) = sym(zeros(BM,1));
for i1=BM+1:TM
    int_rot(i1) = subs(basis(i1),{r,s},{u,v});
    int_rot(i1) = subs(int_rot(i1),{u,v},{s,-1-r-s});
    int_rot2(i1) = subs(int_rot(i1),{r,s},{u,v});
    int_rot2(i1) = subs(int_rot2(i1),{u,v},{s,-1-r-s});
end

% FIND MATRIX COUPLING ROTATED TO NONROTATED
for i1=BM+1:TM
    for i2=BM+1:TM
        int_ms(i1-BM,i2-BM) = int(int(basis(i1)*int_rot(i2),s,-1,-r),r,-1,1);
    end
end

% P = 4 only
% FIND INTERIOR MODES FOR P=4
% syms a2 a3 real
% B = [1-a2-a3,a2,a3]
% t2 = simplify(B*int_ms*B')
% a2 = solve(t2,a2)
% syms a3val;
% a3val = 0.0;
% a2 = subs(a2(1),a3,a3val);
% B = [1-a2-a3val,a2,a3val];
% lumped(BM+1) = B(1)*basis(BM+1) +B(2)*basis(BM+2) +B(3)*basis(BM+3);
% lumped(BM+2) = B(1)*int_rot(BM+1) +B(2)*int_rot(BM+2) +B(3)*int_rot(BM+3);
% lumped(BM+3) = B(1)*int_rot2(BM+1) +B(2)*int_rot2(BM+2) +B(3)*int_rot2(BM+3);
% mode(1) = factor(subs(lumped(BM+1),{r},{(x+1)/2*(1-s)-1})/((x+1)*(x-1)*(s-1)^2*(1+s)))
% mode(2) = factor(subs(lumped(BM+2),{r},{(x+1)/2*(1-s)-1})/((x+1)*(x-1)*(s-1)^2*(1+s)))
% mode(3) = factor(subs(lumped(BM+3),{r},{(x+1)/2*(1-s)-1})/((x+1)*(x-1)*(s-1)^2*(1+s)))
% mode1a = subs(mode,{x,s},{2*u-1,2*v-1})
% pretty(factor(mode1a))

% bdry = (1-x)/2*(1+x)/2*((1-n)/2)^2*(1+n)/2
% mode1 = bdry*(1+ (-10-sqrt(sym(10)))/3*(1+n)/2)
% lumped(BM+1) = simplify(subs(mode1,{x,n},{(2*(1+r)/(1-s) -1),s}));
% v1 = subs(mode1,{x,n},{2*(1+r)/(1-s)-1,s});
% v2 = subs(v1,{r,s},{v,-1-u-v});
% v3 = subs(v2,{u,v},{(1-n)*(1+x)/2-1,n})
% lumped(BM+2) = simplify(subs(v3,{x,n},{(2*(1+r)/(1-s) -1),s}));
% v1 = subs(mode1,{x,n},{2*(1+r)/(1-s)-1,s});
% v2 = subs(v1,{r,s},{-1-u-v,u});
% v3 = subs(v2,{u,v},{(1-n)*(1+x)/2-1,n})
% lumped(BM+3) = simplify(subs(v3,{x,n},{(2*(1+r)/(1-s) -1),s}));


% % TRYING TO FIND INTERIOR MODES FOR P=5
% % SEARCHING FOR SEPARABLE FUNCTIONS...
% % FIND MASS MATRIX OF S FUNCTIONS
% for i=BM+1:TM
%     for j=BM+1:TM
%         sms(i-BM,j-BM) = int(basisn(i)*basisn(j),-1,1);
%     end
% end
% a13 = sms(1:3,1:3)\sms(1:3,6);
% a23 = sms(4:5,4:5)\sms(4:5,6);
% a12 = sms(1:3,1:3)\sms(1:3,4:5);
% 
% % POSSIBLE FUNCTIONS 
% syms a1 a2 a3 a4 a5 a6 b1 b2 b3 b4 b5 b6 real
% B = [1-a2-a3-a4-a5-a6,a2,a3,a4,a5,a6;1-b2-b3-b4-b5-b6,b2,b3,b4,b5,b6]
% B = [a1,a2,a3,a4,a5,a6;b1,b2,b3,b4,b5,b6]
% B = [1-a2-a3,a2,a3,0,0,0;1-b2-b3,b2,b3,0,0,0]; 
% % B = [0,0,0,a4,a5,0;b1,b2,b3,b4,b5,b6];
% % B = [0,0,0,0,0,a6;b1,b2,b3,b4,b5,b6];
% % B = [a1*a13(1),a1*a13(2),a1*a13(3),a2*a23(1),a2*a23(2),1;b1*a13(1),b1*a13(2),b1*a13(3),b2*a23(1),b2*a23(2),1];
% %B = [a1*a12(1,1)+a2*a12(1,2),a1*a12(2,1)+a2*a12(2,2),a1*a12(3,1)+a2*a12(3,2),a3,a4,0;b1*a12(1,1)+b2*a12(1,2),b1*a12(2,1)+b2*a12(2,2),b1*a12(3,1)+b2*a12(3,2),b3,b4,0];
% %B = [a1*a12(1,1)+a2*a12(1,2),a1*a12(2,1)+a2*a12(2,2),a1*a12(3,1)+a2*a12(3,2),a3,a4,0;b1*a13(1),b1*a13(2),b1*a13(3),b2*a23(1),b2*a23(2),b3];
% 
% 
% t1 = B*B';
% t2 = B*int_ms*B';
% eqi(1,1) = t1(2,1);
% eqi(2,1) = t2(1,1);
% eqi(3,1) = t2(2,2);
% eqi(4,1) = t2(2,1);
% eqi(5,1) = t2(1,2);
% eqi = simplify(eqi)
% solve(eqi(2),a3)
% solve(eqi(3),b3)

% NEW WAY TO FIND MODES 
C = int_ms + int_ms.';
[V,E] = eig(C);
for i=1:size(V,1)
    V(:,i) = V(:,i)/sqrt(V(:,i).'*V(:,i));
    % NEED TO SORT EIGENVECTORS IN SOME CONSISTENT WAY OR I WILL GO INSANE
    [Ytemp,indx(i)] = max(abs(double(V(:,i))));
end


[Ytemp,isort] = sort(indx);
V = V(:,isort);
V = simplify(V)

% SYMMETRIC CONSTRAINT
SYMM = E(isort,isort)

% ANTISYMMETRIC CONSTRAINT
ASYM = simplify(V'*(int_ms-int_ms')*V)

% THE FOLLOWING PROVES THERE IS NO SOLUTION USING JUST MODES 1-3
%func = (b1+a3)^2-(2+12/11*a3^2)*(2*b1^2+12/11)
%ezsurf(func,[-3,3,-1,1])
% ALSO NOT POSSIBLE TO USE MODE 6 BY ITSELF


% a1 = input('give me first vector');
% b1 = input('and second');
% a = V*a1;
% b = V*b1;
% 
% double(a'*b)
% double(a'*int_ms*b)
% double(a'*int_ms*a)
% double(b'*int_ms*a)
% double(b'*int_ms*b)
% 
% lumped(BM+1) = a'*basis(BM+1:TM)';
% lumped(BM+2) = b'*basis(BM+1:TM)';
% lumped(BM+3) = subs(lumped(BM+1),{r,s},{u,v});
% lumped(BM+3) = subs(lumped(BM+3),{u,v},{s,-1-r-s});
% lumped(BM+4) = subs(lumped(BM+2),{r,s},{u,v});
% lumped(BM+4) = subs(lumped(BM+4),{u,v},{s,-1-r-s});
% lumped(BM+5) = subs(lumped(BM+3),{r,s},{u,v});
% lumped(BM+5) = subs(lumped(BM+5),{u,v},{s,-1-r-s});
% lumped(BM+6) = subs(lumped(BM+4),{r,s},{u,v});
% lumped(BM+6) = subs(lumped(BM+6),{u,v},{s,-1-r-s});
% 
% % syms u v
% % int_rot(1:BM) = sym(zeros(BM,1));
% for i1=BM+1:TM
%     int_rot(i1) = subs(lumped(i1),{r,s},{u,v});
%     int_rot(i1) = subs(int_rot(i1),{u,v},{s,-1-r-s});
%     int_rot2(i1) = subs(int_rot(i1),{r,s},{u,v});
%     int_rot2(i1) = subs(int_rot2(i1),{u,v},{s,-1-r-s});
% end
% 
% % FIND MATRIX COUPLING ROTATED TO NONROTATED
% for i1=BM+1:TM
%     for i2=BM+1:TM
%         int_ms(i1-BM,i2-BM) = int(int(lumped(i1)*int_rot(i2),s,-1,-r),r,-1,1);
%     end
% end

% LETS SEE ROTATIONAL EIGENFUNCTIONS
lumped(BM+1:TM) = V'*basis(BM+1:TM)';


%for i=BM+1:(TM-BM)/3+mod(TM-BM,3) +BM
for i=BM+1:TM
	figure;
	axis equal;
	zplotm = double(subs(lumped(i),{r,s},{rplotm,splotm}));
    box off;
	contourf(rplotm,splotm,zplotm);
end 

basis = lumped;
matrix_tri;
dms = double(ms);
dms(BM+1:TM,BM+1:TM)

% % ROTATION IN TERMS OF xi AND eta
% v0 = [x,n];
% v1 = subs(v0,{x,n},{2*(1+r)/(1-s)-1,s});
% v2 = subs(v1,{r,s},{v,-1-u-v});
% v3 = subs(v2,{u,v},{(1-n)*(1+x)/2-1,n})
% coord1 = v3;
% 
% v0 = [x,n];
% v1 = subs(v0,{x,n},{2*(1+r)/(1-s)-1,s});
% v2 = subs(v1,{r,s},{-1-u-v,u});
% v3 = subs(v2,{u,v},{(1-n)*(1+x)/2-1,n})
% coord2 = v3;

% syms a b
% figure
% ezcontour(subs(coord1(1),{x,n},{a,b}),[-1,1,-1,1]);
% hold on
% ezcontour(subs(coord1(2),{x,n},{a,b}),[-1,1,-1,1]);
% 
% syms a b
% figure
% ezcontour(subs(coord1(1)+coord2(2),{x,n},{a,b}),[-1,1,-1,1]);
% hold on
% ezcontour(subs(coord1(2)+coord2(1),{x,n},{a,b}),[-1,1,-1,1]);

% % FIND NUMBER OF POSSIBLE UNIQUE SYMMETRIC FUNCTIONS */
% phis(1:BM) = sym(0);
% for i=BM+1:TM
%     phis(i) = basis(i)+int_rot(i)+int_rot2(i);
% end
% 
% % CREATE MASS MATRIX
% for i1=BM+1:TM
%     for i2=BM+1:TM
%         symmas(i1,i2) = int(int(basis(i1)*phis(i2),s,-1,-r),r,-1,1);
%     end
% end
% 
% rank(symmas)
    

        