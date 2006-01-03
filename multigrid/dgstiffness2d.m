Nel = 5;

dx = 0.1/2.0;
dy = 1.0/2.0;

disp(['p = ' num2str(P)]);
cd bases;
clear basis restrict;
switch (basis_flag)
    case 0
        dubiner;
    case 1
        legendre;
    case 2
        lumped;
    case 3
        monomial;
    case 4
        nodal;
    case 5
        intlegendre;
end
matrix2d;

% THIS IS TO TEST MATRIX2D 
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


% THIS IS THE HAROLD TRI/QUAD CASE
% monomial_tri;
% matrix2d_analytic;

cd ..

nvar = 1;
    
if (krestrict & lvl > 1)
    % CAN CALCULATE STIFFNESS MATRIX THIS WAY:
	K = vrestrict{lvl-1}*K*vrestrict{lvl-1}';
	KL = vrestrict{lvl-1}*KL*vrestrict{lvl-1}';
	KR = vrestrict{lvl-1}*KR*vrestrict{lvl-1}';
	KLL = vrestrict{lvl-1}*KLL*vrestrict{lvl-1}';
	KRR = vrestrict{lvl-1}*KRR*vrestrict{lvl-1}';
    KD = vrestrict{lvl-1}*KD*vrestrict{lvl-1}';
	KU = vrestrict{lvl-1}*KU*vrestrict{lvl-1}';
	KDD = vrestrict{lvl-1}*KDD*vrestrict{lvl-1}';
	KUU = vrestrict{lvl-1}*KUU*vrestrict{lvl-1}';
else
    basis2dsave = basis2d;
	if (sys_flag == 0) 
    	%%%%%%%%%%%%%%
		% CONVECTION
		%%%%%%%%%%%%%%
		angle = 10.0*pi/180.0;
		ax = cos(angle);
		ay = sin(angle);
        
        ax = ax*dy;
        ay = ay*dx;
        
        % DETERMINE TAU FOR SUPG UPWINDING
        tau = inv(abs(ax) +abs(ay))/(P+1)^2;
		
		% STIFFNESS MATRICES
		K = ax*cvx +ay*cvy ...
            +(abs(ax)+ax)/2*rl2d -(-abs(ax)+ax)/2*lr2d ...
            +(abs(ay)+ay)/2*tb2d -(-abs(ay)+ay)/2*bt2d;
        
        if (ksupg)  
            K = K + ax*tau*ax*dfx +ay*tau*ay*dfy +ax*tau*ay*dxy +ay*tau*ax*dyx;
        end
		KL = -(abs(ax)+ax)/2*ll2d;
		KR = (-abs(ax)+ax)/2*rr2d;
		KU = (-abs(ay)+ay)/2*tt2d;
		KD = -(abs(ay)+ay)/2*bb2d;
		KLL = zeros(size(KL));
		KRR = zeros(size(KR));
		KDD = zeros(size(KD));
		KUU = zeros(size(KU));
        
	elseif (sys_flag == 1)
		%%%%%%%%%%%%%
		%DIFFUSION %
		%%%%%%%%%%%%%
        if (scheme == 1)
            tstring = 'Central Difference (Bassi-Rebay [10])';
            fstring = 'cd';
		elseif (scheme == 2)
            tstring = 'Brezzi et al. [22]';
            fstring = 'brezzi22';
		elseif (scheme == 3)
            tstring = 'LDG with Central Difference (\beta = 0) [41])';
            fstring = 'ldg_cd';
		elseif (scheme == 4)
            tstring = 'LDG with Upwind Difference (\beta = 1/2) [41])';
            fstring = 'ldg_upwind';
		elseif (scheme == 5)
            tstring = 'Interior Penalty [50]';
            fstring = 'ip50';
		elseif (scheme == 6)
            tstring = 'Bassi et al. [13]';
            fstring = 'bassi13';
		elseif (scheme == 7)
            tstring = 'Upwind Difference';
            fstring = 'upwind';
		end
		
		% CENTRAL DIFFERENCE
		CMX=  (cvx +0.5*rl2d -0.5*lr2d)*dy;
		CL =  -0.5*ll2d*dy;
		CR =  0.5*rr2d*dy;
        CMY = (cvy +0.5*tb2d -0.5*bt2d)*dx;
        CB =  -0.5*bb2d*dx;
		CT =  0.5*tt2d*dx;
		% UPWIND CORRECTION
		UMX = 0.5*(rl2d +lr2d)*dy;
		UL = -0.5*ll2d*dy;
		UR = -0.5*rr2d*dy;
        UMY = 0.5*(tb2d +bt2d)*dx;
        UB = -0.5*bb2d*dx;
		UT = -0.5*tt2d*dx;
        
		if (scheme == 1 || scheme == 2|| scheme == 3)
            % CENTRAL DIFFERENCE
			QMX =  CMX;
			QL =  CL;
			QR =  CR;
            QMY = CMY;
            QB =  CB;
            QT =  CT;
			TMX =  -CMX;
			TL =  -CL;
			TR =  -CR;
            TMY = -CMY;
            TB =  -CB;
            TT =  -CT;
		elseif (scheme == 7 || scheme == 4)
			% UPWINDED
			QMX =  CMX +UMX;
			QL =  CL +UL;
			QR =  CR +UR;
            QMY =  CMY +UMY;
            QB =  CB +UB;
			QT =  CT +UT;
			TMX =  -CMX +UMX;
			TL =  -CL +UL;
			TR =  -CR +UR;
            TMY =  -CMY +UMY;
            TB =  -CB +UB;
			TT =  -CT +UT;
            
		elseif (scheme == 5 || scheme == 6)
			% INTERIOR PENALTY & Bassi [13]
			QMX =  CMX;
			QL =  CL;
			QR =  CR;
            QMY = CMY;
            QB =  CB;
            QT =  CT;
			TMX =  -cvx*dy;
			TL =  zeros(size(TMX));
			TR =  zeros(size(TMX));
            TMY =  -cvy*dx;
            TB =  zeros(size(TMY));
			TT =  zeros(size(TMY));
		end
		
		% PENALTY TERMS
        if (scheme == 3 || scheme == 4) 
			% PENALTY COEFFICENT FOR LDG 
			alphax = 4.0/(2*dx);
            alphay = 4.0/(2*dy);
        elseif (scheme == 5)
            % PENALTY COEFFICIENT FOR IP MUST BE LARGER
            alphax = 20.0/(2*dx);
            alphay = 20.0/(2*dy);
		elseif (scheme == 2 || scheme == 6)
			% PENALTY COEFFICIENT FOR Brezzi et al. [22] / Bassi [13]
			alphax = (((inv(ms)*bright)'*bright +(inv(ms)*bleft)'*bleft)/4)/dx;
            alphay = (((inv(ms)*bright)'*bright +(inv(ms)*bleft)'*bleft)/4)/dy;
            % TO GET P=0 EXACTLY RIGHT
            % alphax = 2*alphax;
            % alphay = 2*alphay;
		else
			alphax = 0;
            alphay = 0;
        end
%        alphax = alphax*4;
%        alphay = alphay*4;
        
		PM = +alphax*(lr2d +rl2d)*dy +alphay*(bt2d +tb2d)*dx;
		PL = -alphax*ll2d*dy;
		PR = -alphax*rr2d*dy;
		PB = -alphay*bb2d*dx;
		PT = -alphay*tt2d*dx;
        
		if (scheme == 5 || scheme == 6)
			% ADD DERIVATIVE TERMS FOR IP / Bassi [13]
			PM = PM +(-1/2*brdr +1/2*bldl)*dy/dx +(-1/2*btdt +1/2*bbdb)*dx/dy;
			PL = PL +1/2*bldr*dy/dx;
			PR = PR -1/2*brdl*dy/dx;
            PB = PB +1/2*bbdt*dx/dy;
			PT = PT -1/2*btdb*dx/dy;
		end
		
		% STATIC INVERSION OF Q
        QMX = inv(dx*dy*ms2d)*QMX;
		QL = inv(dx*dy*ms2d)*QL;
		QR = inv(dx*dy*ms2d)*QR;
        QMY = inv(dx*dy*ms2d)*QMY;
        QB = inv(dx*dy*ms2d)*QB;
        QT = inv(dx*dy*ms2d)*QT;
    
        K = TL*QR +TMX*QMX +TR*QL +TB*QT +TMY*QMY +TT*QB +PM;
		KLL = TL*QL; 
		KL = TL*QMX +TMX*QL +PL;
		KR = TMX*QR +TR*QMX +PR;
		KRR = TR*QR;
        
        KDD = TB*QB; 
		KD = TB*QMY +TMY*QB +PB;
		KU = TMY*QT +TT*QMY +PT;
		KUU = TT*QT;        
	else
        nvar = 4;
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
        
        Ae = Ae*dy;
        absAe = absAe*dy;
        
        Af = Af*dx;
        absAf = absAf*dy;
        
        % DETERMINE TAU FOR SUPG UPWINDING
        tau = inv(absAe +absAf)/(P+1)^2;
		
		% STIFFNESS MATRICES
		K = blktimes(Ae,cvx) +blktimes(Af,cvy) ...
            +blktimes((absAe+Ae)/2,rl2d) -blktimes((-absAe+Ae)/2,lr2d) ...
            +blktimes((absAf+Af)/2,tb2d) -blktimes((-absAf+Af)/2,bt2d);
        
        if (ksupg)  
            K = K + blktimes(Ae*tau*Ae,dfx) +blktimes(Af*tau*Af,dfy) +blktimes(Ae*tau*Af,dxy) +blktimes(Af*tau*Ae,dyx);
        end
		KL = -blktimes((absAe+Ae)/2,ll2d);
		KR =  blktimes((-absAe+Ae)/2,rr2d);
		KU =  blktimes((-absAf+Af)/2,tt2d);
		KD = -blktimes((absAf+Af)/2,bb2d);
		KLL = zeros(size(KL));
		KRR = zeros(size(KR));
		KDD = zeros(size(KD));
		KUU = zeros(size(KU));
   
	end
end

% PRECONDITIONERS
if (sys_flag == 0) 
    switch(rlx_flag)
        case 0
            % MASS MATRIX
	        M = ms2d;
            if (msupg)
                M = M -ax*tau*cvx -ay*tau*cvy;
            end
        case 1
			% JACOBI
			M = diag(diag(K));
        case 2
			% ELEMENT JACOBI
			M = K;
        case 3
            % JACOBI OF inv(ms)*K
            M = inv(inv(diag(diag(inv(ms2d)*K)))*inv(ms2d));
    end
    
    switch(sweep_flag)
        case 0
            % JACOBI TYPE
            ML = zeros(size(M));
            MLL = zeros(size(M));
            MD = zeros(size(M));
            MDD = zeros(size(M));
            MR = zeros(size(M));
            MRR = zeros(size(M));
            MT = zeros(size(M));
            MTT = zeros(size(M));
        case 1
            % GAUSS SEIDEL TYPE
            ML = KL;
            MLL = KLL;
            MD = KD;
            MDD = KDD;
            MR = zeros(size(M));
            MRR = zeros(size(M));
            MT = zeros(size(M));
            MTT = zeros(size(M));
         case 2
            % SYMMETRIC GAUSS SEIDEL TYPE
            ML = KL;
            MLL = KLL;
            MD = KD;
            MDD = KDD;
            MR = KR;
            MRR = KRR;
            MT = KU;
            MTT = KUU;
        case 3
            % LINE SOLVER
            ML = KL;
            MLL = KLL;
            MD = zeros(size(M));
            MDD = zeros(size(M)); 
            MR = KR;
            MRR = KRR;
            MT = zeros(size(M));
            MTT = zeros(size(M));
        case 4
            % SEQUENTIAL LINES
            ML = KL;
            MLL = KLL;
            MD = KD;
            MDD = KDD;
            MR = KR;
            MRR = KRR;
            MT = zeros(size(M));
            MTT = zeros(size(M));
         case 5
            % SYMMETRIC SEQUENTIAL LINES
            ML = KL;
            MLL = KLL;
            MD = KD;
            MDD = KDD;
            MR = KR;
            MRR = KRR;
            MT = zeros(size(M));
            MTT = zeros(size(M));
    end
elseif(sys_flag == 1)
    switch(rlx_flag)
        case 0
            % MASS MATRIX
	        M = ms2d*dx*dy;
        case 1
			% JACOBI
			M = diag(diag(K));
        case 2
			% ELEMENT JACOBI
			M = K;
        case 3
            % JACOBI OF inv(ms)*K
            M = ms2d*dx*dy;
            M = inv( inv( diag( diag(inv(M)*K)) )*inv(dx*dy*ms2d) );
    end
    
    switch(sweep_flag)
        case 0
            % JACOBI TYPE
            ML = zeros(size(M));
            MLL = zeros(size(M));
            MD = zeros(size(M));
            MDD = zeros(size(M));
            MR = zeros(size(M));
            MRR = zeros(size(M));
            MT = zeros(size(M));
            MTT = zeros(size(M));
        case 1
            % GAUSS SEIDEL TYPE
            ML = KL;
            MLL = KLL;
            MD = KD;
            MDD = KDD;
            MR = zeros(size(M));
            MRR = zeros(size(M));
            MT = zeros(size(M));
            MTT = zeros(size(M));
         case 2
            % SYMMETRIC GAUSS SEIDEL TYPE
            ML = KL;
            MLL = KLL;
            MD = KD;
            MDD = KDD;
            MR = KR;
            MRR = KRR;
            MT = KU;
            MTT = KUU;
        case 3
            % LINE SOLVER
            ML = KL;
            MLL = KLL;
            MD = zeros(size(M));
            MDD = zeros(size(M)); 
            MR = KR;
            MRR = KRR;
            MT = zeros(size(M));
            MTT = zeros(size(M));
        case 4
            % SEQUENTIAL LINES
            ML = KL;
            MLL = KLL;
            MD = KD;
            MDD = KDD;
            MR = KR;
            MRR = KRR;
            MT = zeros(size(M));
            MTT = zeros(size(M));
         case 5
            % SYMMETRIC SEQUENTIAL LINES
            ML = KL;
            MLL = KLL;
            MD = KD;
            MDD = KDD;
            MR = KR;
            MRR = KRR;
            MT = zeros(size(M));
            MTT = zeros(size(M));
    end
else
    switch(rlx_flag)
        case 0
            % MASS MATRIX
	        M = blktimes(eye(size(Ae)),ms2d);
            if (msupg)
                M = M -blktimes(Ae*tau,cvx) -blktimes(Af*tau*Af,cvy);
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
    end
    
    switch(sweep_flag)
        case 0
            % JACOBI TYPE
            ML = zeros(size(M));
            MLL = zeros(size(M));
            MD = zeros(size(M));
            MDD = zeros(size(M));
            MR = zeros(size(M));
            MRR = zeros(size(M));
            MT = zeros(size(M));
            MTT = zeros(size(M));
        case 1
            % GAUSS SEIDEL TYPE
            ML = KL;
            MLL = KLL;
            MD = KD;
            MDD = KDD;
            MR = zeros(size(M));
            MRR = zeros(size(M));
            MT = zeros(size(M));
            MTT = zeros(size(M));
         case 2
            % SYMMETRIC GAUSS SEIDEL TYPE
            ML = KL;
            MLL = KLL;
            MD = KD;
            MDD = KDD;
            MR = KR;
            MRR = KRR;
            MT = KU;
            MTT = KUU;
        case 3
            % LINE SOLVER
            ML = KL;
            MLL = KLL;
            MD = zeros(size(M));
            MDD = zeros(size(M)); 
            MR = KR;
            MRR = KRR;
            MT = zeros(size(M));
            MTT = zeros(size(M));
        case 4
            % SEQUENTIAL LINES
            ML = KL;
            MLL = KLL;
            MD = KD;
            MDD = KDD;
            MR = KR;
            MRR = KRR;
            MT = zeros(size(M));
            MTT = zeros(size(M));
         case 5
            % SYMMETRIC SEQUENTIAL LINES
            ML = KL;
            MLL = KLL;
            MD = KD;
            MDD = KDD;
            MR = KR;
            MRR = KRR;
            MT = zeros(size(M));
            MTT = zeros(size(M));
    end
end

% CONSTRUCT A SMALL SQUARE BOX
bsz = size(K,1);
kbig = zeros(Nel*Nel*bsz,Nel*Nel*bsz);
mbig = zeros(Nel*Nel*bsz,Nel*Nel*bsz);

for elx = 1:Nel
    llelx = (elx-3);
    if (llelx < 0) 
        llelx = llelx +Nel;
    end
    lelx = (elx-2);
    if (lelx < 0) 
        lelx = lelx +Nel;
    end
    relx = (elx);
    if (relx > Nel-1) 
        relx = relx-Nel;
    end
    rrelx = (elx+1);
    if (rrelx > Nel-1) 
        rrelx = rrelx-Nel;
    end
    llind = 1 +llelx*bsz:(llelx+1)*bsz;
    lind = 1 +lelx*bsz:(lelx+1)*bsz;
    mindx = 1+(elx-1)*bsz:elx*bsz;
    rind = 1 +(relx)*bsz:(relx+1)*bsz;
    rrind = 1 +(rrelx)*bsz:(rrelx+1)*bsz;
    
    
    for ely = 1:Nel
        llely = (ely-3);
        if (llely < 0) 
            llely = llely +Nel;
        end
        lely = (ely-2);
        if (lely < 0) 
            lely = lely +Nel;
        end
        rely = (ely);
        if (rely > Nel-1) 
            rely = rely-Nel;
        end
        rrely = (ely+1);
        if (rrely > Nel-1) 
            rrely = rrely-Nel;
        end
		llindy = llely*bsz*Nel +mindx;
		lindy = +lely*bsz*Nel +mindx;
		mind = +(ely-1)*bsz*Nel +mindx;
		rindy = +(rely)*bsz*Nel +mindx;
		rrindy = +(rrely)*bsz*Nel +mindx;
        
        llindx = llind +(ely-1)*bsz*Nel;
        lindx = lind +(ely-1)*bsz*Nel;
        rindx = rind +(ely-1)*bsz*Nel;
        rrindx = rrind +(ely-1)*bsz*Nel;

        kbig(mind,mind) = kbig(mind,mind) +K;
        kbig(mind,lindx) = kbig(mind,lindx) +KL;
        kbig(mind,rindx) = kbig(mind,rindx) +KR;
        kbig(mind,llindx) = kbig(mind,llindx) +KLL;
        kbig(mind,rrindx) = kbig(mind,rrindx) +KRR;
        kbig(mind,lindy) = kbig(mind,lindy) +KD;
        kbig(mind,rindy) = kbig(mind,rindy) +KU;
        kbig(mind,llindy) = kbig(mind,llindy) +KDD;
        kbig(mind,rrindy) = kbig(mind,rrindy) +KUU;
        
        mbig(mind,mind) = mbig(mind,mind) +M;
        mbig(mind,lindx) = mbig(mind,lindx) +ML;
        mbig(mind,llindx) = mbig(mind,llindx) +MLL;
        mbig(mind,llindy) = mbig(mind,llindy) +MDD;
        mbig(mind,lindy) = mbig(mind,lindy) +MD;
        mbig(mind,rindx) = mbig(mind,rindx) +MR;
        mbig(mind,rrindx) = mbig(mind,rrindx) +MRR;
        mbig(mind,rindy) = mbig(mind,rindy) +MT;
        mbig(mind,rrindy) = mbig(mind,rrindy) +MTT;
    end
end
if (omega > 0)
	%mtemp = blktimes(eye(Nel*Nel),0.5*(M+M'));
	%maxeig = eigs(kbig,mtemp,1);
	biglam = eig(kbig,-mbig);
	maxeig = max(abs(biglam))
	if (isinf(maxeig))
        biglam = sort(biglam);
        count = size(biglam,1);
        while (isinf(biglam(count)))
            count = count -1;
        end
        biglam = biglam(1:count);
        maxeig = biglam(count);
	end
	% figure;
	% plot(real(biglam),imag(biglam),'+');
else
    maxeig = 1.0;
end
maxeig = abs(maxeig)

% STORE ROWS FOR ELEMENT 3 OF STIFFNESS MATRIX FOR FOURIER ANALYSIS
rows = (2*Nel +2)*bsz+1:(2*Nel +3)*bsz;
for m=1:5
    for n=1:5
        KB{lvl,n,m} = kbig(rows,rows +(m-3)*bsz +(n-3)*Nel*bsz);
        MB{lvl,n,m} = mbig(rows,rows +(m-3)*bsz +(n-3)*Nel*bsz);
    end
end
vmaxeig(lvl) = maxeig;
vrestrict{lvl} = blktimes(eye(nvar,nvar),restrict2d);

% FOR DIFFUSION AT P=1, BLOCK JACOBI, WITH NO SWEEPING
if (P == 1  && sys_flag == 1  && rlx_flag == 2 && (sweep_flag == 0 || sweep_flag == 3)) 
    omega = 2/3*omega
end

% RESCALE DIAGONAL TERM OF PRECONDITIONER
if (sweep_flag < 3)
    % FOR NON-LINE SOLVES
    MB{lvl,3,3} = MB{lvl,3,3}/abs(omega)*maxeig;
else
    % FOR LINE SOLVES
    for m = 1:5
        MB{lvl,3,m} = MB{lvl,3,m}/abs(omega)*maxeig;
    end
end
    

for m=1:3
    for n=1:3
        RB{lvl,n,m} = zeros(nvar*size(restrict2d));
    end
end
RB{lvl,2,2} = blktimes(eye(nvar,nvar),restrict2d);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TESTING TO CHECK RESULTS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if (lvl == 1)
%    Nel = 8;
% 
%     % CONSTRUCT STIFFNESS MATRIX
%     Kbigx = zeros(Nel*size(KB{lvl,3,3},1),Nel*size(KB{lvl,3,3},1));
%     Mbigx = zeros(size(Kbigx));
%     Rbigx = zeros(Nel*size(RB{lvl,2,2},1),Nel*size(RB{lvl,2,2},2));
%     for m = 1:5
%         position = (m-2);
%         if (position <= 0) 
%             position = position +Nel;
%         end
%         pvect = zeros(1,Nel);
%         pvect(position) = 1;
%         Kbigx = Kbigx + blktimes(gallery('circul',pvect),KB{lvl,3,m});
%         Mbigx = Mbigx + blktimes(gallery('circul',pvect),MB{lvl,3,m});
%     end
%     for m = 1:3
%         position = (m-1);
%         if (position <= 0) 
%             position = position +Nel;
%         end
%         pvect = zeros(1,Nel);
%         pvect(position) = 1;
%         Rbigx = Rbigx + blktimes(gallery('circul',pvect),RB{lvl,2,m});
%     end
% 
%     % PLACE X MATRICES ON DIAGONAL
%     pvect = zeros(1,Nel);
%     pvect(1) = 1;
%     Kbig = blktimes(gallery('circul',pvect),Kbigx);
%     Mbig = blktimes(gallery('circul',pvect),Mbigx);
%     Rbig = blktimes(gallery('circul',pvect),Rbigx);
% 
%     for m = 1:5
%         if (m==3)
%             continue
%         end
%         pvect = zeros(1,Nel);
%         pvect(1) = 1;
%         Kbigy = blktimes(gallery('circul',pvect),KB{lvl,m,3});
%         Mbigy = blktimes(gallery('circul',pvect),MB{lvl,m,3});
%         position = (m-2);
%         if (position <= 0) 
%             position = position +Nel;
%         end
%         pvect = zeros(1,Nel);
%         pvect(position) = 1;
%         Kbig = Kbig + blktimes(gallery('circul',pvect),Kbigy);
%         Mbig = Mbig + blktimes(gallery('circul',pvect),Mbigy);
%     end    
% 
% 
%     % TEST TO SEE IF I'VE GOT THIS RIGHT
%     Kbigc = Rbig*Kbig*(Rbig');
%     fine = (eye(size(Kbig)) -inv(Mbig)*Kbig);
%     fineeig = eig(fine);
%     wholething = (eye(size(Kbig)) -Rbig'*lscov(Kbigc,Rbig*Kbig))*fine;
%     
%     eigswt =eig(wholething);
%     figure
%     plot(real(fineeig),imag(fineeig),'x');
%     title('amplification factor (direct)');
%     axis equal
%     afactordirect =max(abs(fineeig(2:size(fineeig,1))))
% 
%     figure
%     plot(real(eigswt),imag(eigswt),'x')
%     title('multigrid amplification factor (direct)');
%     axis equal
%     dfactordirect =max(abs(eigswt(2:size(eigswt,1))))
% end


if (~continuous_flag) 
    return;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOLLOWING LINES ARE FOR COARSENING TO CONTINUOUS SPACE AFTER P = 1 %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('switching to continuous space');
    
bsz = size(K,1);
nvar = bsz/(P+1)^2;
bszc = nvar*P^2;
rbig = zeros(size(kbig,1)/4,size(kbig,2));

% COARSEN TO CONTINUOUS SPACE
% ASSUMES MODES ARE ORGANIZED IN STRIPS IN X (NOT CCW)
for elx = 1:Nel
    lelx = elx-1;
    if (lelx == 0)
        lelx = Nel;
    end
    lind = 1 +(lelx-1)*bsz:lelx*bsz;
    mindx = 1+(elx-1)*bsz:elx*bsz;
    mindxc = 1+(elx-1)*bszc:elx*bszc;
    
    for ely = 1:Nel
        lely = ely-1;
        if (lely == 0)
            lely = Nel;
        end
		lindy = (lely-1)*bsz*Nel +mindx;
		mind = (ely-1)*bsz*Nel +mindx;
        lindx = lind +(ely-1)*bsz*Nel;
        lindxy = lind +(lely-1)*bsz*Nel;
        
		mindc = +(ely-1)*bszc*Nel +mindxc;
        
        for indx=1:nvar
            rbig(mindc,lindx+1:(P+1)^2:lindx+bsz-1) = 1.0;
            rbig(mindc,mind:(P+1)^2:mind+bsz-1) = 1.0;
            rbig(mindc,lindy+2:(P+1)^2:lindy+bsz-1) = 1.0;
            rbig(mindc,lindxy+3:(P+1)^2:lindxy+bsz-1) = 1.0;
        end
    end
end


% HAROLD'S IDEA OF TRYING TO FORCE IT TO BE CELL-CENTERED
% rbig2 = zeros(size(kbig,1)/4,size(kbig,1)/4);
% for elx = 1:Nel
%     lelx = elx-1;
%     if (lelx == 0)
%         lelx = Nel;
%     end
%     lind = 1 +(lelx-1)*bszc:lelx*bszc;
%     mindx = 1+(elx-1)*bszc:elx*bszc;
%     
%     for ely = 1:Nel
%         lely = ely-1;
%         if (lely == 0)
%             lely = Nel;
%         end
% 		lindy = (lely-1)*bszc*Nel +mindx;
%         lindxy = (lely-1)*bszc*Nel +lind;
% 
% 		mind = mindx +(ely-1)*bszc*Nel;
%         lindx = lind +(ely-1)*bszc*Nel;
%         
% 		mindc = +(ely-1)*bszc*Nel +mindx;
%         
%         for indx=1:nvar
%             rbig2(mind,lindx) = 1.0;
%             rbig2(mind,mind) = 1.0;
%             rbig2(mind,lindy) = 1.0;
%             rbig2(mind,lindxy) = 1.0;
%         end
%     end
% end
% rbig = rbig2*rbig;

% FORM COARSE GRID STIFFNESS MATRIX
kbigc = rbig*kbig*rbig.';

% STORE ROWS FOR ELEMENT 3 OF STIFFNESS MATRIX FOR FOURIER ANALYSIS
rows = (2*Nel +2)*bszc+1:(2*Nel +3)*bszc;
for m=1:5
    for n=1:5
        KB{lvl,n,m} = kbigc(rows,rows +(m-3)*bszc +(n-3)*Nel*bszc);
    end
end


rows = (2*Nel +2)*bszc+1:(2*Nel +3)*bszc;
for m=1:3
    for n=1:3
        cols = (n*Nel +m)*bsz+1:(n*Nel +(m+1))*bsz;
        RB{lvl-1,n,m} = rbig(rows,cols);
    end
end


% % TEST TO SEE IF I'VE GOT THIS RIGHT
% fine = (eye(size(kbig)) -inv(mbig)*kbig);
% [finevect,fineeig] = eig(fine);
% wholething = (eye(size(kbig)) -rbig'*lscov(kbigc,rbig*kbig))*fine;
% %wholething = (eye(size(kbig)) -rbig'*inv(kbigc)*rbig*kbig)*fine;
% eigs=eig(wholething);
% figure
% plot(real(fineeig),imag(fineeig),'x');
% title('amplification factor (direct)');
% axis equal
% figure
% plot(real(eigs),imag(eigs),'x')
% title('multigrid amplification factor (direct)');
% axis equal
% dfactordirect=max(abs(eigs(2:size(eigs,1))))