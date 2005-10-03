Nel = 32;
rswp = 0;
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
cd ..

if (krestrict & lvl > 1)
    % CAN CALCULATE STIFFNESS MATRIX THIS WAY:
	K = vrestrict{lvl-1}*K*vrestrict{lvl-1}';
	KL = vrestrict{lvl-1}*KL*vrestrict{lvl-1}';
	KR = vrestrict{lvl-1}*KR*vrestrict{lvl-1}';
	KLL = vrestrict{lvl-1}*KLL*vrestrict{lvl-1}';
	KRR = vrestrict{lvl-1}*KRR*vrestrict{lvl-1}';
    
else
    basissave = basis;
    
	%%%%%%%%%%%%%%%
	% CONVECTION:
	%%%%%%%%%%%%%%
	if (sys_flag == 0) 
		% STIFFNESS MATRICES
		K = cv +rightl;
		KL = -leftl;
		KR = zeros(size(KL));
		% ADD SUPG UPWINDING
		if (ksupg) 
            K = K +1.0/(P+1)^2*df;
        end
		KLL = zeros(size(K));
		KRR = zeros(size(K));
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
		% elseif (scheme == 9)
		%     tstring = 'Implict Boundary Flux';
		%     fstring = 'iflx';
		end
		if (rswp)
            tstring = [tstring ' with Swapping'];
            fstring = [fstring 'swp'];
		end
		
		% CENTRAL DIFFERENCE
		CM =  cv +0.5*rightl -0.5*leftr;
		CL =  -0.5*leftl;
		CR =  0.5*rightr;
		% UPWIND CORRECTION
		UM = 0.5*(rightl +leftr);
		UL = -0.5*leftl;
		UR = -0.5*rightr;
		
		if (scheme == 1 || scheme == 2|| scheme == 3)
            % CENTRAL DIFFERENCE
			QM =  CM;
			QL =  CL;
			QR =  CR;
			TM =  -CM;
			TL =  -CL;
			TR =  -CR;
		elseif (scheme == 7 || scheme == 4)
			% UPWINDED
			QM =  CM +UM;
			QL =  CL +UL;
			QR =  CR +UR;
			TM =  -CM +UM;
			TL =  -CL +UL;
			TR =  -CR +UR;
		elseif (scheme == 5 || scheme == 6)
			% INTERIOR PENALTY & Bassi [13]
			QM =  CM;
			QL =  CL;
			QR =  CR;
			TM =  -cv;
			TL =  zeros(size(cv));
			TR =  zeros(size(cv));
		% elseif (scheme == 9)
		%     % IMPLICT BOUNDARY HEAT FLUX
		%     QM = cv +rightl -leftr;
		%     QL = zeros(size(cv));
		%     QR = zeros(size(cv));
		%     TM = -cv;
		%     TL = zeros(size(cv));
		%     TR = zeros(size(cv));
		end
		
		% % NO STATIC INVERSION OF Q
		% K = [ms,QM;TM,zeros(size(TM))];
		% KL = [zeros(size(QM)),QL;TL,zeros(size(TM))];
		% KR = [zeros(size(QM)),QR;TR,zeros(size(TM))];
		% K = zeros(size(K));
		% KRR = zeros(size(KR));
		% % PRECONDITIONER
		% % JUST MASS MATRIX
		% M = [ms,zeros(size(ms));zeros(size(ms)),ms];
		% restrict = [restrict,zeros(size(restrict));zeros(size(restrict)),restrict];
		
		% PENALTY TERMS
		if (scheme == 3 || scheme == 4) 
			% PENALTY COEFFICENT FOR LDG 
			alpha = 4.0/2;
        elseif (scheme == 5)
            % PENALTY COEFFICIENT FOR IP MUST BE LARGER
            alpha = 20.0/2;
		elseif (scheme == 2 || scheme == 6)
			% PENALTY COEFFICIENT FOR Brezzi et al. [22] / Bassi [13]
			alpha = ((inv(ms)*bright)'*bright +(inv(ms)*bleft)'*bleft)/4;
            alpha = 2*alpha;  % TO GET P=0 EXACTLY RIGHT
		else
			alpha = 0;
		end
        
        %alpha = alpha/2;
        
		PM = +alpha*leftr +alpha*rightl;
		PL = -alpha*leftl;
		PR = -alpha*rightr;
        
        alpha
		
		if (scheme == 5 || scheme == 6)
			% ADD DERIVATIVE TERMS FOR IP / Bassi [13]
			PM = PM -1/2*bright*bdright' +1/2*bleft*bdleft';
			PL = PL +1/2*bleft*bdright';
			PR = PR -1/2*bright*bdleft';
		end
		
		% STATIC INVERSION OF Q
		QL = inv(ms)*QL;
		QM = inv(ms)*QM;
		QR = inv(ms)*QR;
		KLL = TL*QL; 
		KL = TL*QM +TM*QL +PL;
		K = TL*QR +TM*QM +TR*QL +PM;
		KR = TM*QR +TR*QM +PR;
		KRR = TR*QR;
        
        % PREMULTIPLY BY inv(M)?
%         KLL = inv(ms)*KLL;
%         KL = inv(ms)*KL;
%         K = inv(ms)*K;
%         KR = inv(ms)*KR;
%         KRR = inv(ms)*KRR;
		    
		if (rswp)
			for iedge = 1:Nel
                % lrswitch(iedge) = 2*floor(2*rand(1,1))-1;
                % lrswitch(iedge) = 1;
                lrswitch(iedge) = 2*floor(2*(iedge-1)/Nel)-1;
			end
			lrswitch(Nel+1) = lrswitch(1);
            for iedge = 1:Nel
                lrswitch(iedge) = 0.5*(lrswitch(iedge) + lrswitch(iedge+1));
            end
            lrswitch(Nel+1) = lrswitch(1);
            
            lrswitch = lrswitch*0.5;
		end
	else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% LINEARIZED EULER (quasi-1D)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        K =  blktimes(Ae,cv) +blktimes((absAe+Ae)/2,rightl) -blktimes((-absAe+Ae)/2,leftr);
		if (ksupg) 
            K = K +blktimes(absAe,1.0/(P+1)^2*df);
        end
        KL = -blktimes((absAe+Ae)/2,leftl);
		KR =  blktimes((-absAe+Ae)/2,rightr);
		KLL = zeros(size(K));
		KRR = zeros(size(K));
	
        restrict = blktimes(eye(size(Ae)),restrict);
	end
end

if (sys_flag == 0)
% PRECONDITIONERS
	% MASS MATRIX
    switch(rlx_flag)
        case 0
			% ADD SUPG UPWINDING?
			if (msupg) 
                M = M -1.0/(P+1)^2*cv;
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
            if (msupg) 
                M = M - 2.0/(P+1)^2*cv;
            end
            M = inv(inv(diag(diag(inv(M)*K)))*inv(M));
    end
    
    switch(sweep_flag)
        case 0
            % JACOBI TYPE
            ML = zeros(size(M));
            MLL = zeros(size(M));
            MR = zeros(size(M));
            MRR = zeros(size(M));
        case 1
            % GAUSS SEIDEL TYPE
            ML = KL;
            MLL = KLL;
            MR = zeros(size(M));
            MRR = zeros(size(M));
        case 2
            % SYMMETRIC GAUSS SEIDEL TYPE
            ML = KL;
            MLL = KLL;
            %M = M +KL*inv(M)*KR;
            MR = KR;
            MRR = KRR;
    end
elseif (sys_flag == 1)
    % PRECONDITIONER
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
            M = inv(inv(diag(diag(inv(ms)*K)))*inv(ms));

    end
    
    switch(sweep_flag)
        case 0
            % JACOBI TYPE
            ML = zeros(size(M));
            MLL = zeros(size(M));
            MR = zeros(size(M));
            MRR = zeros(size(M));
        case 1
            % GAUSS SEIDEL TYPE
            ML = KL;
            MLL = KLL;
            MR = zeros(size(M));
            MRR = zeros(size(M));
        case 2
            % SYMMETRIC GAUSS SEIDEL TYPE
            ML = KL;
            MLL = KLL;
            %M = M +KL*inv(M)*KR;
            M = M;
            MR = KR;
            MRR = KRR;
    end
elseif (sys_flag == 2)
    % PRECONDITIONERS
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
            if (msupg) 
                M = M -blktimes(Ve*sign(Ee)*inv(Ve),2.0/(P+1)^2*cv);
            end
            M = inv(inv(diag(diag(inv(M)*K)))*inv(M));
    end
    
    switch(sweep_flag)
        case 0
            % JACOBI TYPE
            ML = zeros(size(M));
            MLL = zeros(size(M));
            MR = zeros(size(M));
            MRR = zeros(size(M));
        case 1
            % GAUSS SEIDEL TYPE
            ML = KL;
            MLL = KLL;
            MR = zeros(size(M));
            MRR = zeros(size(M));
        case 2
            % SYMMETRIC GAUSS SEIDEL TYPE
            ML = KL;
            MLL = KLL;
            MR = KR;
            MRR = KRR;
    end
end

% CALCULATE MAGNITUDE OF EIGENVALUES FOR RELAXATION SCHEME
bsz = size(K,1);
kbig = zeros(Nel*bsz,Nel*bsz);
mbig = zeros(Nel*bsz,Nel*bsz);
for el=1:Nel
    llel = (el-3);
    if (llel < 0) 
        llel = llel +Nel;
    end
    lel = (el-2);
    if (lel < 0) 
        lel = lel +Nel;
    end
    rel = (el);
    if (rel > Nel-1) 
        rel = rel-Nel;
    end
    rrel = (el+1);
    if (rrel > Nel-1) 
        rrel = rrel-Nel;
    end
    llind = 1 +llel*bsz:(llel+1)*bsz;
    lind = 1 +lel*bsz:(lel+1)*bsz;
    mind = 1+(el-1)*bsz:el*bsz;
    rind = 1 +rel*bsz:(rel+1)*bsz;
    rrind = 1 +rrel*bsz:(rrel+1)*bsz;
    
    if (rswp) 
        % SWAPPING TEST
        QM =  CM +lrswitch(el)*leftr +lrswitch(el+1)*rightl;
		QL =  CL -lrswitch(el)*leftl;
		QR =  CR -lrswitch(el+1)*rightr;
		TM =  -CM +lrswitch(el)*leftr +lrswitch(el+1)*rightl;
		TL =  -CL -lrswitch(el)*leftl;
		TR =  -CR -lrswitch(el+1)*rightr;
        QM = inv(ms)*QM;
        QL = inv(ms)*QL;
        QR = inv(ms)*QR;
		KLL = TL*QL; 
		KL = TL*QM +TM*QL +PL;
		K = TL*QR +TM*QM +TR*QL +PM;
		KR = TM*QR +TR*QM +PR;
		KRR = TR*QR;
        
        % IF GAUSS SEIDEL & SWAPPING
        % ML = KL;
        % MLL = KLL;
    end
    
    kbig(mind,mind) = kbig(mind,mind) +K;
    kbig(mind,lind) = kbig(mind,lind) +KL;
    kbig(mind,rind) = kbig(mind,rind) +KR;
    kbig(mind,llind) = kbig(mind,llind) +KLL;
    kbig(mind,rrind) = kbig(mind,rrind) +KRR;
    
    mbig(mind,mind) = mbig(mind,mind) +M;
    if (omega < 0)
       mbig(mind,lind) = mbig(mind,lind) +ML;
       mbig(mind,llind) = mbig(mind,llind) +MLL;
       mbig(mind,rind) = mbig(mind,rind) +MR;
       mbig(mind,rrind) = mbig(mind,rrind) +MRR;
    end
end

if (omega < 0) 
    maxeig = 1;
else
    biglam = eig(kbig,mbig);
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
end

%figure;
%plot(real(biglam),imag(biglam),'+');

% STORE ROWS FOR ELEMENT 3 OF STIFFNESS MATRIX FOR FOURIER ANALYSIS
rows = 2*bsz+1:3*bsz;
for m=1:5
    KB{lvl,m} = kbig(rows,rows+(m-3)*bsz);
end
MB{lvl,1} = MLL;
MB{lvl,2} = ML;
MB{lvl,3} = M;
MB{lvl,4} = MR;
MB{lvl,5} = MRR;

% RESCALE PRECONDITIONER DIAGONAL BY MAXEIG
MB{lvl,3} = MB{lvl,3}/abs(omega)*maxeig;
vmaxeig{lvl} = maxeig;
vrestrict{lvl} = restrict;

for m=1:3
    RB{lvl,m} = zeros(size(restrict));
end
RB{lvl,2} = restrict;

if (rswp) 
    %figure;
    %spy(kbig);
    rvec = sort(-real(biglam));
    figure;
    subplot(2,1,1)
    plot(rvec,'x');
    wnmber(1) = 0;
    for k = 2:size(rvec,1)
        wnmber(k) = 2.*pi/ (2*Nel/(floor(k/2)));
        rvec(k) = rvec(k)/wnmber(k)^2;
    end
    subplot(2,1,2)
    plot(wnmber(2:size(wnmber,2)),rvec(2:size(rvec,1)),'x');
    return;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TESTING TO CHECK RESULTS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCT STIFFNESS MATRIX
% if (lvl == 1)
%     Nel = 7;
%     Kbig = zeros(Nel*size(KB{lvl,3},1),Nel*size(KB{lvl,3},1));
%     Mbig = zeros(size(Kbig));
%     Rbig = zeros(Nel*size(RB{lvl,2},1),Nel*size(RB{lvl,2},2));
%     for m = 1:5
%         position = (m-2);
%         if (position <= 0) 
%             position = position +Nel;
%         end
%         pvect = zeros(1,Nel);
%         pvect(position) = 1;
%         Kbig = Kbig + blktimes(gallery('circul',pvect),KB{lvl,m});
%         Mbig = Mbig + blktimes(gallery('circul',pvect),MB{lvl,m});
%     end
%     for m = 1:3
%         position = (m-1);
%         if (position <= 0) 
%             position = position +Nel;
%         end
%         pvect = zeros(1,Nel);
%         pvect(position) = 1;
%         Rbig = Rbig + blktimes(gallery('circul',pvect),RB{lvl,m});
%     end
%     % TEST TO SEE IF I'VE GOT THIS RIGHT
%     Kbigc = Rbig*Kbig*(Rbig')
%     fine = (eye(size(Kbig)) -inv(Mbig)*Kbig);
%     fineeig = eig(fine);
%     rfactordirect =max(abs(fineeig))
% 
%     wholething = (eye(size(Kbig)) -Rbig'*lscov(Kbigc,Rbig*Kbig))*fine;
%     eigswt =eig(wholething);
%     figure
%     plot(real(fineeig),imag(fineeig),'x');
%     title('amplification factor (direct)');
%     axis equal
%     figure
%     plot(real(eigswt),imag(eigswt),'x')
%     title('multigrid amplification factor (direct)');
%     axis equal
%     dfactordirect =max(abs(eigswt(2:size(eigswt,1))))
% end




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOLLOWING LINES ARE FOR COARSENING TO CONTINUOUS SPACE AFTER P = 1 %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~continuous_flag) 
    return;
end


disp('switching to continuous space');
    
nvar = size(K,1)/(P+1);
bsz = nvar*P;
rbig = zeros(size(kbig,1)/2,size(kbig,2));

% COARSEN TO CONTINUOUS SPACE
for el=1:Nel
    elm = el-1;
    if (elm == 0) 
        elm = Nel;
    end
    for indx=1:nvar
        rbig((el-1)*nvar+indx,2*(el-1)*nvar+(indx-1)*2+1) = 1.0;
        rbig((el-1)*nvar+indx,2*(elm-1)*nvar+(indx-1)*2+2) = 1.0;
    end
end
kbigc = rbig*kbig*rbig.';
kbigc
%mbigc = rbig*mbig*rbig.';
%biglam = eig(-inv(mbig)*kbig);
% figure;
% plot(real(biglam),imag(biglam),'+');
% maxeig = max(abs(biglam))

% STORE ROWS FOR ELEMENT 3 OF STIFFNESS MATRIX FOR FOURIER ANALYSIS
% THIS ONLY WORKS FOR NODAL BASES AND SCHEMES WITH BANDWITH 3
rows = 2*bsz+1:3*bsz;
for m=1:5
    KB{lvl,m} = kbigc(rows,rows+(m-3)*bsz);
%    MB{lvl,m} = mbigc(rows,rows+(m-3)*bsz);
end
%vmaxeig{lvl} = maxeig;
%vrestrict{lvl} = restrict;

RB{lvl-1,1} = rbig(nvar+1:2*nvar,1:2*nvar);
RB{lvl-1,2} = rbig(nvar+1:2*nvar,2*nvar+1:4*nvar);
RB{lvl-1,3} = rbig(nvar+1:2*nvar,4*nvar+1:6*nvar);


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
% dfactor=max(abs(eigs(2:size(eigs,1))))



