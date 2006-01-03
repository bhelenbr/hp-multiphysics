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
matrix2d;
cd ..

if (lvl > 1 & krestrict)
    % CAN CALCULATE STIFFNESS MATRIX THIS WAY:
	K = vrestrict{lvl-1}*K*vrestrict{lvl-1}';
else
	%%§%%%%%%%%%%%%
	% CONVECTION:
	%%%%%%%%%%%%%%
	if (sys_flag == 0) 
        angle = 0.0*pi/180.0;
		ax = cos(angle);
		ay = sin(angle);
		% STIFFNESS MATRICES
		% FOR SUPG UPWINDING 
		K = ax*cvx +ay*cvy;
		if (ksupg)
            tau = 2/(2*1*(P+1)^2);
			K = K +tau*(ax^2*dfx +ay^2*dfy +ax*ay*dxy +ay*ax*dyx);
		end
	elseif (sys_flag == 1)
	%%%%%%%%%%%%%
	%DIFFUSION %
	%%%%%%%%%%%%%	
		K = dfx +dfy;
	else
        %%%%%%%%%%%%%%%%%%%
        % LINEARIZED EULER
        %%%%%%%%%%%%%%%%%%%
		syms w1 w2 w3 w4 positive;
		syms gam positive;
		e = [w2;w2^2/w1+(gam-1)*(w4-1/2*(w2^2+w3^2)/w1);w2*w3/w1;w2/w1*(w4+(gam-1)*(w4-1/2*(w2^2+w3^2)/w1))];
		f = [w3;w2*w3/w1;w3^2/w1+(gam-1)*(w4-1/2*(w2^2+w3^2)/w1);w3/w1*(w4+(gam-1)*(w4-1/2*(w2^2+w3^2)/w1))];
        Ae = jacobian(e,[w1,w2,w3,w4]);
        Af = jacobian(f,[w1,w2,w3,w4]);
	
        syms rho u v gam c positive
        Ae = subs(Ae,{w1,w2,w3,w4},{rho,rho*u,rho*v,rho*(c^2/(gam*(gam-1)) +1/2*(u^2+v^2))});
        Af = subs(Af,{w1,w2,w3,w4},{rho,rho*u,rho*v,rho*(c^2/(gam*(gam-1)) +1/2*(u^2+v^2))});
        [Ve,Ee] = eig(Ae);
        absAe = simplify(Ve*abs(Ee)*inv(Ve));
		[Vf,Ef] = eig(Af);
        absAf = simplify(Vf*abs(Ef)*inv(Vf));
	
        mach = 0.1;
        gam = 1.4;
        angle = 0.0*pi/180.0;
        rho = 1;
        u = mach*cos(angle);
        v = mach*sin(angle);
        c = 1;
		Ae = double(subs(Ae));
        absAe = double(subs(absAe));
        Af = double(subs(Af));
		absAf = double(subs(absAf));
        disp('finished in symbolic land');
		
		% STIFFNESS MATRICES
		K = blktimes(Ae,cvx) +blktimes(Af,cvy)
		if (ksupg) 
            tau = 2/(2*(P+1)^2)*inv(absAe +absAf);
            K = K +blktimes(Ae*tau*Ae,dfx) +blktimes(Af*tau*Af)*dfy ...
                  +blktimes(Ae*tau*Af,dxy) +blktimes(Af*tau*Ae,dyx);
        end
        
        restrict = blktimes(eye(size(Ae)),restrict);
	end
end

% PRECONDITIONERS
if (sys_flag == 0)
    switch(rlx_flag)
        case 0
            % MASS MATRIX
	        M = ms2d;
			% ADD SUPG UPWINDING?
			if (msupg)
				M = M - ax*tau*cvx -ay*tau*cvy;
			end
        case 1
			% JACOBI
			M = diag(diag(K));
        case 2
			% ELEMENT JACOBI
			M = K;
        case 3
            % JACOBI OF inv(ms)*K
             M = ms2d;
			% ADD SUPG UPWINDING?
			if (msupg)
				M = M -ax*tau*cvx -ay*tau*cvy;
			end
            M = inv(inv(diag(diag(inv(M)*K)))*inv(M));
        case 4
        	M = mapp2d;
        case 5
        	M = inv(inv(diag(diag(inv(mapp2d)*K)))*inv(mapp2d));
    end
elseif(sys_flag == 1)
    switch(rlx_flag)
        case 0
	        M = ms2d;
        case 1
			% JACOBI
			M = diag(diag(K));
        case 2
			% ELEMENT JACOBI
			M = K;
        case 3
            % JACOBI OF inv(ms)*K
            M = ms2d;
            M = inv(inv(diag(diag(inv(M)*K)))*inv(M));
        case 4
        	M = mapp2d;
        case 5
        	M = inv(inv(diag(diag(inv(mapp2d)*K)))*inv(mapp2d));
    end
else
    switch(rlx_flag)
        case 0
			M = blktimes(eye(size(Ae)),ms2d);
			% ADD SUPG UPWINDING?
			if (msupg) 
                M = M -blktimes(Ae*tau,cvx) -blktimes(Af*tau*cvy);
            end
        case 1
			% JACOBI
			M = diag(diag(K));
        case 2
			% ELEMENT JACOBI
			M = K;
        case 3
            % JACOBI OF inv(ms)*K
            M = blktimes(eye(size(Ae)),ms2d);
            M = inv(inv(diag(diag(inv(M)*K)))*inv(M));
        case 4
        	M = blktimes(eye(size(Ae)),mapp2d);
        case 5
            M = blktimes(eye(size(Ae)),mapp2d);
        	M = inv(inv(diag(diag(inv(M)*K)))*inv(M));
    end
end

% BRUTE FORCE WAY AS CHECK
nvar = size(K,1)/(P+1)^2;
bsz = nvar*P^2;
kbig = zeros(Nel*Nel*bsz,Nel*Nel*bsz);
mbig = zeros(Nel*Nel*bsz,Nel*Nel*bsz);

% SET UP INDEXING
% INTERIOR
indx = [];
for m = 1:P
    rowb = (m-1)*P;
    indx = [indx rowb+1:rowb+P rowb+1+bsz];
end
indx = [indx Nel*bsz+1:Nel*bsz+P Nel*bsz+bsz+1];

% RIGHT SIDE
indxr = [];
base = (Nel-1)*bsz;
for m = 1:P
    rowb = base +(m-1)*P;
    indxr = [indxr rowb+1:rowb+P (m-1)*P+1];
end
indxr = [indxr base+Nel*bsz+1:base+Nel*bsz+P Nel*bsz+1];

% TOP SIDE
indxt = [];
base = (Nel-1)*Nel*bsz;
for m = 1:P
    rowb = base +(m-1)*P;
    indxt = [indxt rowb+1:rowb+P rowb+1+bsz];
end
indxt = [indxt 1:P bsz+1];

% UPPER RIGHT ELEMENT
indxur = [];
base = ((Nel-1)*Nel +(Nel-1))*bsz;
for m = 1:P
    rowb = base +(m-1)*P;
    indxur = [indxur rowb+1:rowb+P (Nel-1)*Nel*bsz+1+(m-1)*P];
end
indxur = [indxur (Nel-1)*bsz+1:(Nel-1)*bsz+P 1];

for ely = 1:Nel-1
    for elx = 1:Nel-1
		eind = +((ely-1)*Nel +(elx-1))*bsz;
        for n=1:nvar
            rows = indx+eind+(n-1)*P^2;
            for m=1:nvar
                cols = indx+eind+(m-1)*P^2;
                kbig(rows,cols) = kbig(rows,cols) +K((n-1)*(P+1)^2+1:n*(P+1)^2,(m-1)*(P+1)^2+1:m*(P+1)^2);
                mbig(rows,cols) = mbig(rows,cols) +M((n-1)*(P+1)^2+1:n*(P+1)^2,(m-1)*(P+1)^2+1:m*(P+1)^2);
            end
        end
    end
    % RIGHT SIDE
    eind = (ely-1)*Nel*bsz;
    for n=1:nvar
        rows = indxr+eind+(n-1)*P^2;
        for m=1:nvar
            cols = indxr+eind+(m-1)*P^2;
            kbig(rows,cols) = kbig(rows,cols) +K((n-1)*(P+1)^2+1:n*(P+1)^2,(m-1)*(P+1)^2+1:m*(P+1)^2);
            mbig(rows,cols) = mbig(rows,cols) +M((n-1)*(P+1)^2+1:n*(P+1)^2,(m-1)*(P+1)^2+1:m*(P+1)^2);
        end
    end    
end
% TOP SIDE
for elx = 1:Nel-1
    eind = (elx-1)*bsz;
    for n=1:nvar
        rows = indxt+eind+(n-1)*P^2;
        for m=1:nvar
            cols = indxt+eind+(m-1)*P^2;
            kbig(rows,cols) = kbig(rows,cols) +K((n-1)*(P+1)^2+1:n*(P+1)^2,(m-1)*(P+1)^2+1:m*(P+1)^2);
            mbig(rows,cols) = mbig(rows,cols) +M((n-1)*(P+1)^2+1:n*(P+1)^2,(m-1)*(P+1)^2+1:m*(P+1)^2);
        end
    end
end
% UPPER RIGHT ELEMENT
for n=1:nvar
    rows = indxur+(n-1)*P^2;
    for m=1:nvar
        cols = indxur+(m-1)*P^2;
        kbig(rows,cols) = kbig(rows,cols) +K((n-1)*(P+1)^2+1:n*(P+1)^2,(m-1)*(P+1)^2+1:m*(P+1)^2);
        mbig(rows,cols) = mbig(rows,cols) +M((n-1)*(P+1)^2+1:n*(P+1)^2,(m-1)*(P+1)^2+1:m*(P+1)^2);
    end
end
biglam = eig(-inv(mbig)*kbig);
% figure;
% plot(real(biglam),imag(biglam),'+');
maxeig = max(abs(biglam))

% STORE ROWS FOR ELEMENT 3 OF STIFFNESS MATRIX FOR FOURIER ANALYSIS
rows = (2*Nel +2)*bsz+1:(2*Nel +3)*bsz;
for m=1:5
    for n=1:5
        KB{lvl,n,m} = kbig(rows,rows +(m-3)*bsz +(n-3)*Nel*bsz);
        MB{lvl,n,m} = mbig(rows,rows +(m-3)*bsz +(n-3)*Nel*bsz);
    end
end
vmaxeig(lvl) = maxeig;

if (lvl < maxlvl)
	% ASSEMBLE RESTRICTION OPERATOR 
	% FINE MESH
	indx = [];
	for m = 1:P
        rowb = (m-1)*P;
        indx = [indx rowb+1:rowb+P rowb+1+bsz];
	end
	indx = [indx Nel*bsz+1:Nel*bsz+P Nel*bsz+bsz+1];
	
	% COARSE MESH
	indxc = [];
	PC = sqrt(size(restrict2d,1)/nvar) -1;
	bszc = nvar*PC^2;
	for m = 1:PC;
        rowb = (m-1)*PC;
        indxc = [indxc rowb+1:rowb+PC rowb+1+bszc];
	end
	indxc = [indxc Nel*bszc+1:Nel*bszc+PC Nel*bszc+bszc+1];
	
	% NOT COMPLETE DIDN'T DO ENDPOINTS
	rbig = zeros(Nel*Nel*bszc,Nel*Nel*bsz);
	for ely = 1:Nel-1
        for elx = 1:Nel-1
			eind = +((ely-1)*Nel +(elx-1))*bsz;
            eindc = ((ely-1)*Nel +(elx-1))*bszc;
            for n=1:nvar
                rows = indxc+eindc+(n-1)*PC^2;
                for m=1:nvar
                    cols = indx+eind+(m-1)*P^2;
                    rbig(rows,cols) = restrict2d((n-1)*(PC+1)^2+1:n*(PC+1)^2,(m-1)*(P+1)^2+1:m*(P+1)^2);
                end
            end
        end
	end
	rows = (2*Nel +2)*bsz+1:(2*Nel +3)*bsz;
	crows = (2*Nel +2)*bszc+1:(2*Nel +3)*bszc;
	vrestrict{lvl} = restrict2d;
    for m=1:3
        for n=1:3
            RB{lvl,n,m} = rbig(crows,rows +(m-2)*bsz +(n-2)*Nel*bsz);
        end
	end
end

% % UNCOMMENT SPECIAL CASE OF GEOMETRIC MULTGRID 
% TKB(:,:) = KB(1,:,:);
% for m=1:5
%     for n=1:5
%         KB{1,m,n} = zeros(4,4);
%     end
% end
% KB{1,3,3} = [TKB{3,3},TKB{3,4},TKB{4,3},TKB{4,4};TKB{3,2},TKB{3,3},TKB{4,2},TKB{4,3};TKB{2,3},TKB{2,4},TKB{3,3},TKB{3,4};TKB{2,2},TKB{2,3},TKB{3,2},TKB{3,3}];
% KB{1,3,2} = [0,TKB{3,2},0,TKB{4,2};0,0,0,0;0,TKB{2,2},0,TKB{3,2};0,0,0,0];
% KB{1,3,4} = [0,0,0,0;TKB{3,4},0,TKB{4,4},0;0,0,0,0;TKB{2,4},0,TKB{3,4},0];
% KB{1,2,2} = [0,0,0,TKB{2,2};0,0,0,0;0,0,0,0;0,0,0,0];
% KB{1,2,3} = [0,0,TKB{2,3},TKB{2,4};0,0,TKB{2,2},TKB{2,3};0,0,0,0;0,0,0,0];
% KB{1,2,4} = [0,0,0,0;0,0,TKB{2,4},0;0,0,0,0;0,0,0,0];
% KB{1,4,2} = [0,0,0,0;0,0,0,0;0,TKB{4,2},0,0;0,0,0,0];
% KB{1,4,3} = [0,0,0,0;0,0,0,0;TKB{4,3},TKB{4,4},0,0;TKB{4,2},TKB{4,3},0,0];
% KB{1,4,4} = [0,0,0,0;0,0,0,0;0,0,0,0;TKB{4,4},0,0,0];
% 
% % FOR DIFFUSION
% for m=1:5
%     for n=1:5
%         KB{2,m,n} = TKB{m,n};
%     end
% end
% 
% TKB(:,:) = MB(1,:,:);
% for m=1:5
%     for n=1:5
%         MB{1,m,n} = zeros(4,4);
%     end
% end
% MB{1,3,3} = [TKB{3,3},TKB{3,4},TKB{4,3},TKB{4,4};TKB{3,2},TKB{3,3},TKB{4,2},TKB{4,3};TKB{2,3},TKB{2,4},TKB{3,3},TKB{3,4};TKB{2,2},TKB{2,3},TKB{3,2},TKB{3,3}];
% MB{1,3,2} = [0,TKB{3,2},0,TKB{4,2};0,0,0,0;0,TKB{2,2},0,TKB{3,2};0,0,0,0];
% MB{1,3,4} = [0,0,0,0;TKB{3,4},0,TKB{4,4},0;0,0,0,0;TKB{2,4},0,TKB{3,4},0];
% MB{1,2,2} = [0,0,0,TKB{2,2};0,0,0,0;0,0,0,0;0,0,0,0];
% MB{1,2,3} = [0,0,TKB{2,3},TKB{2,4};0,0,TKB{2,2},TKB{2,3};0,0,0,0;0,0,0,0];
% MB{1,2,4} = [0,0,0,0;0,0,TKB{2,4},0;0,0,0,0;0,0,0,0];
% MB{1,4,2} = [0,0,0,0;0,0,0,0;0,TKB{4,2},0,0;0,0,0,0];
% MB{1,4,3} = [0,0,0,0;0,0,0,0;TKB{4,3},TKB{4,4},0,0;TKB{4,2},TKB{4,3},0,0];
% MB{1,4,4} = [0,0,0,0;0,0,0,0;0,0,0,0;TKB{4,4},0,0,0];
% 
% RB{1,1,1} = [0,0,0,0.25];
% RB{1,2,1} = [0,0.5,0.0,0.25];
% RB{1,3,1} = [0,0,0,0];
% RB{1,1,2} = [0,0,.5,0.25];
% RB{1,2,2} = [1,0.5,0.5,0.25];
% RB{1,3,2} = [0,0,0,0];
% RB{1,1,3} = [0,0,0,0];
% RB{1,2,3} = [0,0,0,0];
% RB{1,3,3} = [0,0,0,0];