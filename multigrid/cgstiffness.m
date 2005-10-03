Nel = 5;

P = max(P,1);
PC = max(floor(P/2),1);

disp(['p = ' num2str(P)]);
cd bases;
clear basis restrict;
switch (basis_flag)
    case 0
        dubiner;
    case 1
        disp('can not use legendre');
    case 2
        lumped;
    case 3
        disp('can not use monomial');
    case 4
        nodal;
end
cd ..

if (lvl > 1 & krestrict)
    % CAN CALCULATE STIFFNESS MATRIX THIS WAY:
	K = vrestrict{lvl-1}*K*vrestrict{lvl-1}';
else
	%%§%%%%%%%%%%%%
	% CONVECTION:
	%%%%%%%%%%%%%%
	if (sys_flag == 0) 
		% STIFFNESS MATRICES
		% FOR SUPG UPWINDING 
		K = cv;
		if (ksupg)
			K = K +1.0/(P+1)^2*df;
		end
	elseif (sys_flag == 1)
	%%%%%%%%%%%%%
	%DIFFUSION %
	%%%%%%%%%%%%%	
		K = df;
	else
		% LINEARIZED EULER (quasi-1D)
		syms w1 w2 w3 w4
		syms gam positive;
		e = [w2;w2^2/w1+(gam-1)*(w4-1/2*(w2^2+w3^2)/w1);w2*w3/w1;w2/w1*(w4+(gam-1)*(w4-1/2*(w2^2+w3^2)/w1))];
		Ae = jacobian(e,[w1,w2,w3,w4]);
		mach = 0.1;
        gam = 1.4;
        w1 = 1;
        w2 = mach;
        w3 = 0.0;
        w4 = 1/(gam*(gam-1)) +1/2*mach^2;
		Ae = subs(Ae);
		Ae = double(Ae);
		[Ve,Ee] = eig(Ae);
		absAe = Ve*abs(Ee)*inv(Ve);
	
		% STIFFNESS MATRICES
        K =  blktimes(Ae,cv);
		if (ksupg) 
            K = K +blktimes(absAe,1.0/(P+1)^2*df);
        end
        
        restrict = blktimes(eye(size(Ae)),restrict);
	end
end

% PRECONDITIONERS
if (sys_flag == 0)
    switch(rlx_flag)
        case 0
            % MASS MATRIX
	        M = ms;
			% ADD SUPG UPWINDING?
			if (msupg)
				M = M - 1.0/(P+1)^2*cv;
			end
        case 1
			% JACOBI
			M = diag(diag(K));
        case 2
			% ELEMENT JACOBI
			M = K;
        case 3
            % JACOBI OF inv(ms)*K
             M = ms;
			% ADD SUPG UPWINDING?
			if (msupg)
				M = M - 1.0/(P+1)^2*cv;
			end
            M = inv(inv(diag(diag(inv(M)*K)))*inv(M));
        case 4
        	M = mapp;
        case 5
        	M = inv(inv(diag(diag(inv(mapp)*K)))*inv(mapp));
    end
elseif(sys_flag == 1)
    switch(rlx_flag)
        case 0
	        M = ms;
        case 1
			% JACOBI
			M = diag(diag(K));
        case 2
			% ELEMENT JACOBI
			M = K;
        case 3
            % JACOBI OF inv(ms)*K
            M = ms;
            M = inv(inv(diag(diag(inv(M)*K)))*inv(M));
        case 4
        	M = mapp;
        case 5
        	M = inv(inv(diag(diag(inv(mapp)*K)))*inv(mapp));
    end
else
    switch(rlx_flag)
        case 0
			M = blktimes(eye(size(Ae)),ms);
			% ADD SUPG UPWINDING?
			if (msupg) 
                M = M -blktimes(Ve*sign(Ee)*inv(Ve),1.0/(P+1)^2*cv);
            end
        case 1
			% JACOBI
			M = diag(diag(K));
        case 2
			% ELEMENT JACOBI
			M = K;
        case 3
            % JACOBI OF inv(ms)*K
            M = blktimes(eye(size(Ae)),ms);
            M = inv(inv(diag(diag(inv(M)*K)))*inv(M));
        case 4
        	M = blktimes(eye(size(Ae)),mapp);
        case 5
            M = blktimes(eye(size(Ae)),mapp);
        	M = inv(inv(diag(diag(inv(M)*K)))*inv(M));
    end
end

% BRUTE FORCE WAY AS CHECK
nvar = size(K,1)/(P+1);
bsz = nvar*P;
kbig = zeros(Nel*bsz,Nel*bsz);
mbig = zeros(Nel*bsz,Nel*bsz);
for el=1:Nel-1
    indx = (el-1)*bsz;
    for n=1:nvar
        rows = [indx+(n-1)*P+1:indx+n*P,(indx+1 +(n-1)*P+bsz)];
        for m=1:nvar
            cols = [indx+(m-1)*P+1:indx+m*P,indx+(m-1)*P+bsz+1];
            kbig(rows,cols) = kbig(rows,cols) +K((n-1)*(P+1)+1:n*(P+1),(m-1)*(P+1)+1:m*(P+1));
            mbig(rows,cols) = mbig(rows,cols) +M((n-1)*(P+1)+1:n*(P+1),(m-1)*(P+1)+1:m*(P+1));
        end
    end
end;
el = Nel;
indx = (el-1)*bsz;
for n=1:nvar
    rows = [(indx+(n-1)*P +1):(indx+n*P),((n-1)*P +1)];
    for m=1:nvar
        cols = [(indx+(m-1)*P +1):(indx+m*P),((m-1)*P +1)];
        kbig(rows,cols) = kbig(rows,cols) +K((n-1)*(P+1)+1:n*(P+1),(m-1)*(P+1)+1:m*(P+1));
        mbig(rows,cols) = mbig(rows,cols) +M((n-1)*(P+1)+1:n*(P+1),(m-1)*(P+1)+1:m*(P+1));
    end
end
biglam = eig(-inv(mbig)*kbig);
% figure;
% plot(real(biglam),imag(biglam),'+');
maxeig = max(abs(biglam))

% STORE ROWS FOR ELEMENT 3 OF STIFFNESS MATRIX FOR FOURIER ANALYSIS
rows = 2*bsz+1:3*bsz;
for m=1:5
    KB{lvl,m} = kbig(rows,rows+(m-3)*bsz);
    MB{lvl,m} = mbig(rows,rows+(m-3)*bsz);
end
vmaxeig(lvl) = maxeig;
vrestrict{lvl} = restrict;

r1 = restrict(1:PC+1,1:P+1);
RB{lvl,1} = blktimes(eye(nvar,nvar),[r1(PC+1,1:P);zeros(PC-1,P)]);
RB{lvl,2} = blktimes(eye(nvar,nvar),r1(1:PC,1:P));
RB{lvl,3} = blktimes(eye(nvar,nvar),[r1(1:PC,P+1),zeros(PC,P-1)]);


% % UNCOMMENT SPECIAL CASE OF GEOMETRIC MULTGRID 
% left = KB{lvl,2};
% mid = KB{lvl,3};
% right = KB{lvl,4};
% KB{1,1} = zeros(2,2);
% KB{1,2} = [0,left;0,0];
% KB{1,3} = [mid,right;left,mid];
% KB{1,4} = [0,0;right,0];
% KB{1,5} = zeros(2,2);
% 
% % FOR DIFFUSION
% KB{2,1} = 0;
% KB{2,2} = left*0.5;
% KB{2,3} = mid*0.5;
% KB{2,4} = right*0.5;
% KB{2,5} = 0;
% 
% left = MB{lvl,2};
% mid = MB{lvl,3};
% right = MB{lvl,4};
% MB{1,1} = zeros(2,2);
% MB{1,2} = [0,left;0,0];
% MB{1,3} = [mid,right;left,mid];
% MB{1,4} = [0,0;right,0];
% MB{1,5} = zeros(2,2);
% 
% RB{1,1} = [0,0.5];
% RB{1,2} = [1,0.5];
% RB{1,3} = [0,0];






