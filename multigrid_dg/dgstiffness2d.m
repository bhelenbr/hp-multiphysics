dx = 1.0/2.0;  % THIS IS Dx/Dxi
dy = dyodx*dx;  % THIS IS Dy/Deta

if (P==0 && lvl==1 && maxlvl > 1)
    dx=dx/2; 
    dy=dy/2;
end

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
if (oned)
    matrix1d;
else
    matrix2d;
end
bsz = size(ms2d,1);

if (lvl == 1)
    basissave = basis2d;
end

% THIS IS THE HAROLD TRI/QUAD CASE
% disp('Harold tri test');
% monomial_tri;
% matrix2d_analytic;
cd ..
    
if (krestrict && lvl > 1) 
    
    % CAN CALCULATE STIFFNESS MATRIX THIS WAY:
% 	KB{lvl,3,3} = RB{lvl-1,2,2}*KB{lvl-1,3,3}*RB{lvl-1,2,2}';
% 	KB{lvl,3,2} = RB{lvl-1,2,2}*KB{lvl-1,3,2}*RB{lvl-1,2,2}';
% 	KB{lvl,3,4} = RB{lvl-1,2,2}*KB{lvl-1,3,4}*RB{lvl-1,2,2}';
%     KB{lvl,3,1} = RB{lvl-1,2,2}*KB{lvl-1,3,1}*RB{lvl-1,2,2}';
%     KB{lvl,3,5} = RB{lvl-1,2,2}*KB{lvl-1,3,5}*RB{lvl-1,2,2}';
%     KB{lvl,2,3} = RB{lvl-1,2,2}*KB{lvl-1,2,3}*RB{lvl-1,2,2}';
% 	KB{lvl,4,3} = RB{lvl-1,2,2}*KB{lvl-1,4,3}*RB{lvl-1,2,2}';
% 	KB{lvl,1,3} = RB{lvl-1,2,2}*KB{lvl-1,1,3}*RB{lvl-1,2,2}';
% 	KB{lvl,5,3} = RB{lvl-1,2,2}*KB{lvl-1,5,3}*RB{lvl-1,2,2}';
%         
    % MATRIX MULTIPLICATION OF RB*KB FOR ELEMENT DISTRIBUTED RB (CONTINUOUS)
    % FOR DG USUALLY ONLY RB(lvl,2,2) is nonzero making this massive overkill
    % IF BANDWITH OF KB IS 5, LARGE BANDWITH TERMS WILL BE TRUNCATED
    for m=1:5
        for n=1:5
            RBKB{m,n} = zeros(size(RB{lvl-1,2,2}*KB{lvl-1,3,3}));
            for rm=1:3
                for rn=1:3
                    offsety = m-(rm-2);
                    offsetx = n-(rn-2);
                    if (offsetx < 1 || offsetx > 5)
                        continue;
                    end
                    if (offsety < 1 || offsety > 5)
                        continue;
                    end
                    RBKB{m,n} = RBKB{m,n} +RB{lvl-1,rm,rn}*KB{lvl-1,offsety,offsetx};
                end
            end
        end
    end
 
    % NOW MATRIX MULTIPLICATION OF (RB*KB)*RB'
    % WHICH IS EQUAL TO (RB*(RB*KB)')' SO CAN USE SAME ALGORITHM AS ABOVE
    % THEN TRANSPOSE
    for m=1:5
        for n=1:5
            RKR{m,n} = zeros(size(RB{lvl-1,2,2}*(RBKB{3,3}')));
            for rm=1:3
                for rn=1:3
                    offsety = 6 -(m-(rm-2));
                    offsetx = 6 -(n-(rn-2));
                    if (offsetx < 1 || offsetx > 5)
                        continue;
                    end
                    if (offsety < 1 || offsety > 5)
                        continue;
                    end
                    RKR{m,n} = RKR{m,n} +RB{lvl-1,rm,rn}*(RBKB{offsety,offsetx}');
                end
            end
        end
    end    
    
    % AND THE FINAL TRANSPOSE
    for m=1:5
        for n=1:5
            offsety = (6-m);
            offsetx = (6-n);
            KB{lvl,m,n} = RKR{offsety,offsetx}';
        end
    end
    
else    
    nvar = 1;
	if (sys_flag == 1 || sys_flag == 2) 
    	%%%%%%%%%%%%%%
		% CONVECTION
		%%%%%%%%%%%%%%
		ax = cos(angle);
		ay = sin(angle);
        if (sys_flag == 2)
            ax = ax*Pe;
            ay = ay*Pe;
        end
        ax = ax*dy;
        ay = ay*dx;
        
        % DETERMINE TAU FOR SUPG UPWINDING
        tau = inv(abs(ax) +abs(ay))/(0.25*(P+1)^2);  % 1/2 isoparametric element length of 2
        if (sys_flag == 2)
            tau = tau*(coth(Pe)-1/Pe); 
        end
        KDX = ax*cvx +(abs(ax)+ax)/2*rl2d -(-abs(ax)+ax)/2*lr2d;
        KDY = +ay*cvy +(abs(ay)+ay)/2*tb2d -(-abs(ay)+ay)/2*bt2d;

		% DIAGONAL
		KB{lvl,3,3} =  KDX +KDY;
        
        if (ksupg)  
            KDX = KDX + ax*tau*ax*dfx;
            KDY = KDY + ay*tau*ay*dfy;
            KB{lvl,3,3} = KB{lvl,3,3} + ax*tau*ax*dfx +ay*tau*ay*dfy +ax*tau*ay*dxy +ay*tau*ax*dyx;
        end
        
        % HORIZONTAL
        KB{lvl,3,1} = zeros(size(KB{lvl,3,3}));
		KB{lvl,3,2} = -(abs(ax)+ax)/2*ll2d;
		KB{lvl,3,4} = (-abs(ax)+ax)/2*rr2d;
        KB{lvl,3,5} = zeros(size(KB{lvl,3,3}));
        

        % VERTICAL
        KB{lvl,1,3} = zeros(size(KB{lvl,3,3}));
		KB{lvl,2,3} = -(abs(ay)+ay)/2*bb2d;
		KB{lvl,4,3} = (-abs(ay)+ay)/2*tt2d;
		KB{lvl,5,3} = zeros(size(KB{lvl,3,3}));
    end
    
    
    
	if (sys_flag == 0 || sys_flag == 2)
		%%%%%%%%%%%%%
		%DIFFUSION %
		%%%%%%%%%%%%%
        if (sys_flag == 0) 
              for m=1:5
                  for n=1:5
                    KB{lvl,m,n} = zeros(bsz,bsz);
                  end
              end
              KDX = zeros(bsz,bsz);
              KDY = zeros(bsz,bsz);
        end
        
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
			alphax = eta*4.0/(2*dx);
            alphay = eta*4.0/(2*dy);
        elseif (scheme == 5)
            % PENALTY COEFFICIENT FOR IP MUST BE LARGER
            alphax = eta*20.0/(2*dx);
            alphay = eta*20.0/(2*dy);
		elseif (scheme == 2 || scheme == 6)
			% PENALTY COEFFICIENT FOR Brezzi et al. [22] / Bassi [13]
			alphax = eta*(((inv(ms)*bright)'*bright +(inv(ms)*bleft)'*bleft)/4)/dx;
            alphay = eta*(((inv(ms)*bright)'*bright +(inv(ms)*bleft)'*bleft)/4)/dy;
            % TO GET P=0 EXACTLY RIGHT
            % alphax = 2*alphax;
            % alphay = 2*alphay;
		else
			alphax = 0;
            alphay = 0;
        end

        eta_actual = alphax*2*dx
        
		PM = +alphax*(lr2d +rl2d)*dy +alphay*(bt2d +tb2d)*dx;
        PMX = +alphax*(lr2d +rl2d)*dy;
        PMY = +alphay*(bt2d +tb2d)*dx;
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
    
        % DIAGONAL
        KDX = KDX +TL*QR +TMX*QMX +TR*QL +PMX;
        KDY = KDY +TB*QT +TMY*QMY +TT*QB +PMY;
        KB{lvl,3,3} = KB{lvl,3,3} +TL*QR +TMX*QMX +TR*QL +TB*QT +TMY*QMY +TT*QB +PM;
		
        % HORIZONTAL
        KB{lvl,3,1} = KB{lvl,3,1} +TL*QL;
		KB{lvl,3,2} = KB{lvl,3,2} +TL*QMX +TMX*QL +PL;
		KB{lvl,3,4} = KB{lvl,3,4} +TMX*QR +TR*QMX +PR;
		KB{lvl,3,5} = KB{lvl,3,5} +TR*QR;
        
        % VERTICAL 
        KB{lvl,1,3} = KB{lvl,1,3} +TB*QB; 
		KB{lvl,2,3} = KB{lvl,2,3} +TB*QMY +TMY*QB +PB;
		KB{lvl,4,3} = KB{lvl,4,3} +TMY*QT +TT*QMY +PT;
		KB{lvl,5,3} = KB{lvl,5,3} +TT*QT;        
    elseif (sys_flag == 3)
        nvar = 4;
        %%%%%%%%%%%%%%%%%%%
        % LINEARIZED EULER
        %%%%%%%%%%%%%%%%%%%
        gam = 1.4;
        rho = 1;
        u = mach*cos(angle);
        v = mach*sin(angle);
        c = 1;
        
        efix = 0.0; % Entropy fix
        if (~prcndtn)
            % UPWINDED STIFFNESS MATRIX USING CONSERVATIVE VARIABLES
            Ee = diag([u-c, u+c,   u,    u]);
            Ve =[1,1,1,0;-c+u,c+u,u,0;v,v,0,1;1/2*(-u^2+2*u*c+u^2*gam+v^2*gam-2*c*u*gam-v^2+2*c^2)/(gam-1),1/2*(u^2*gam-u^2+2*c*u*gam-2*u*c+2*c^2+v^2*gam-v^2)/(gam-1),1/2*(u^2-v^2),v];
            ax = Ve*Ee*inv(Ve);
            absax = Ve*(abs(Ee) +efix*c*mach*eye(4,4))*inv(Ve);

            Ef = diag([v-c, v+c,   v,    v]);
            Vf =[1,1,1,0;u,u,0,1;v-c,v+c,v,0;1/2*(-2*c*v*gam+2*c*v+2*c^2-u^2+v^2*gam-v^2+u^2*gam)/(gam-1),1/2*(2*c*v*gam-2*c*v+2*c^2-u^2+v^2*gam-v^2+u^2*gam)/(gam-1),-1/2*(u^2-v^2),u];
            ay = Vf*Ef*inv(Vf);
            absay = Vf*(abs(Ef) +efix*c*mach*eye(4,4))*inv(Vf);
        else
            % THIS IS THE STIFFNESS MATRIX EXPRESSED USING SYMMETRY VARIABLES
            ax =[u, c, 0, 0; c, u, 0, 0; 0, 0, u, 0; 0, 0, 0, u];
            ay =[v, 0, c, 0; 0, v, 0, 0; c, 0, v, 0; 0, 0, 0, v];
            
            % SIMPLEST PRECONDITIONER 
            % (SEE TURKEL 1999 Ann. Review Fluid Mech.)
            beta2 = (u^2+v^2)/c^2;
            precon = diag([beta2,1,1,1]);
            
            ax = precon*ax;
            speed = sqrt(((1+beta2)*u)^2 +4*beta2*(c^2-u^2));
            lam1 = 1/2*((1+beta2)*u +speed);
            lam2 = 1/2*((1+beta2)*u -speed);
            Ee = diag([lam1,lam2,u,u]);
            Ve = [(beta2-1)*u+speed,(beta2-1)*u-speed,0,0;2*c,2*c,0,0;0,0,1,0;0,0,0,1];
            absax = Ve*(abs(Ee) +efix*c*mach*eye(4,4))*inv(Ve);
    
            ay = precon*ay;
            speed = sqrt(((1+beta2)*v)^2 +4*beta2*(c^2-v^2));
            lam1 = 1/2*((1+beta2)*v +speed);
            lam2 = 1/2*((1+beta2)*v -speed);
            Ef = diag([lam1,lam2,v,v]);
            Vf = [(beta2-1)*v+speed,(beta2-1)*v-speed,0,0;0,0,1,0;2*c,2*c,0,0;0,0,0,1];
            absay = Vf*(abs(Ef) +efix*c*mach*eye(4,4))*inv(Vf);
        end

        ax = ax*dy;
        absax = absax*dy;
        
        ay = ay*dx;
        absay = absay*dy;
        
        if (oned)
            ay = zeros(size(ax));
            absay = zeros(size(ax));
        end
        
        % DETERMINE TAU FOR SUPG UPWINDING
        tau = inv(absax +absay)/(0.25*(P+1)^2);  % 1/2 dpsi = 1 for element length
		
		% DIAGONAL OF STIFFNESS MATRIX    ...LINEARIZED EULER...
		KB{lvl,3,3} = blktimes(ax,cvx) +blktimes(ay,cvy) ...
            +blktimes((absax+ax)/2,rl2d) -blktimes((-absax+ax)/2,lr2d) ...
            +blktimes((absay+ay)/2,tb2d) -blktimes((-absay+ay)/2,bt2d);
        
        KDX = blktimes(ax,cvx) +blktimes((absax+ax)/2,rl2d) -blktimes((-absax+ax)/2,lr2d);
        KDY = blktimes(ay,cvy) +blktimes((absay+ay)/2,tb2d) -blktimes((-absay+ay)/2,bt2d);
        
        if (ksupg)  
            KDX = KDX + blktimes(ax*tau*ax,dfx);
            KDY = KDY + blktimes(ay*tau*ay,dfy);
            KB{lvl,3,3} = KB{lvl,3,3} + blktimes(ax*tau*ax,dfx) +blktimes(ay*tau*ay,dfy) +blktimes(ax*tau*ay,dxy) +blktimes(ay*tau*ax,dyx);
        end
        
        % HORIZONTAL
        KB{lvl,3,1} = zeros(size(KB{lvl,3,3}));
		KB{lvl,3,2} = -blktimes((absax+ax)/2,ll2d);
		KB{lvl,3,4} =  blktimes((-absax+ax)/2,rr2d);
        KB{lvl,3,5} = zeros(size(KB{lvl,3,3}));

        % VERTICAL
        KB{lvl,1,3} = zeros(size(KB{lvl,3,3}));
		KB{lvl,2,3} = -blktimes((absay+ay)/2,bb2d);
		KB{lvl,4,3} =  blktimes((-absay+ay)/2,tt2d);
		KB{lvl,5,3} = zeros(size(KB{lvl,3,3}));
    end

    % ZERO UNUSED BLOCKS
    for m=1:5
        for n=1:5
            if (m ~= 3 && n ~= 3) 
                KB{lvl,m,n} = zeros(size(KB{lvl,3,3}));
            end
        end
    end
end

% DON'T NEED RELAXATION SCHEME ON COARSEST LEVEL
if (lvl == maxlvl && maxlvl > 1) 
    return;
end

% ZERO PRECONDITIONERS BY DEFAULT
for m=1:5
    for n=1:5
        MB{lvl,m,n} = zeros(size(KB{lvl,3,3}));
        MB1{lvl,m,n} = zeros(size(KB{lvl,3,3}));
    end
end

% DIAGONAL BLOCKS OF RELAXATION MATRIX
hblocks = [3,3]; % HORIZONTAL BLOCKS TO USE WHEN CALCULATING MAXEIG
switch(rlx_flag)
    case 0
        % MASS MATRIX
        MB{lvl,3,3} = blktimes(eye(nvar,nvar),ms2d)*dx*dy;
        if (msupg) 
            MB{lvl,3,3} = MB{lvl,3,3} -ax*tau*cvx -ay*tau*cvy;
        end
    case 1
        % JACOBI
        MB{lvl,3,3} = diag(diag(KB{lvl,3,3}));
    case 2
        % ELEMENT JACOBI
        MB{lvl,3,3} = KB{lvl,3,3};
    case 3
        % JACOBI OF M^-1 K
        MB{lvl,3,3} = ms2d*diag(diag(blktimes(eye(nvar,nvar),ms2d)\KB{lvl,3,3}));
    case 4
        % LINE SOLVE
        for k=1:5
            MB{lvl,3,k} = KB{lvl,3,k};
        end
        hblocks = [1,5];
    case 5
        % ALTERNATING DIRECTION LINE SOLVE
        for k=1:5
            MB{lvl,3,k} = KB{lvl,3,k};
            MB1{lvl,k,3} = KB{lvl,k,3};
        end
        hblocks = [1,5];
    case 6
        % ALTERNATING DIRECTION IMPLICIT
        % OFF DIAGONAL TERMS
        for k=[1,2,4,5]
            MB{lvl,3,k} = mu*dt*((dx*dy*blktimes(eye(nvar,nvar),ms2d))\KB{lvl,3,k});
            MB1{lvl,k,3} = mu*dt*((dx*dy*blktimes(eye(nvar,nvar),ms2d))\KB{lvl,k,3});
        end
        
        % DIAGONAL TERMS
        MB{lvl,3,3} = eye(nvar*size(ms2d)) +mu*dt*((dx*dy*blktimes(eye(nvar,nvar),ms2d))\KDX);
        MB1{lvl,3,3} = eye(nvar*size(ms2d)) +mu*dt*((dx*dy*blktimes(eye(nvar,nvar),ms2d))\KDY);
        
        for k=1:5
            MB{lvl,3,k} = dx*dy*blktimes(eye(nvar,nvar),ms2d)*MB{lvl,3,k}/dt;
        end  
    case 7
        % INCOMPLETE L-U
        
        % KTEMP HAS MODIFIED DIAGONAL BLOCK FOR IMPLICIT TIME ADVANCEMENT/UNDERRELAXATION
        for m=2:4
            for n=2:4
                KTEMP{m,n} = KB{lvl,m,n};
            end
        end

        if (implicit_flag)
            for m=2:4
                for n=2:4
                    KTEMP{m,n} = mu*KB{lvl,m,n};
                end
            end
            KTEMP{3,3} = KTEMP{3,3} +blktimes(eye(nvar,nvar),ms2d)/dt;
        end
        
        % 2D PERIODIC
        % ITERATION TO FIND DIAGONAL BLOCK
        if (Fourier2D)

            MB1{lvl,3,3} = KTEMP{3,3};

            Numiter = 1000;
            for n=1:Numiter
                Mnew = KTEMP{3,3} -KTEMP{3,2}*(MB1{lvl,3,3}\KTEMP{3,4}) -KTEMP{2,3}*(MB1{lvl,3,3}\KTEMP{4,3});
                fpterror = norm(Mnew -MB1{lvl,3,3});
				if (fpterror < 1.0e-8) 
					break;
				end
				MB1{lvl,3,3} = Mnew;
            end
            if (n > Numiter-1) 
                disp('uh-oh diagonal not converged in ILU');
            end
            MB{lvl,3,3} = eye(size(MB1{lvl,3,3}));
            MB{lvl,3,2} = KTEMP{3,2}/MB1{lvl,3,3};
            MB{lvl,2,3} = KTEMP{2,3}/MB1{lvl,3,3};

            MB1{lvl,3,4} = KTEMP{3,4};
            MB1{lvl,4,3} = KTEMP{4,3};

            % TEST FOR SCALAR CASE
%             MB1{lvl,3,3}
%             a = KTEMP{2,3}
%             b = KTEMP{3,2}
%             c = KTEMP{3,3}
%             d = KTEMP{3,4}
%             e = KTEMP{4,3}
%             D = (c + sqrt(c^2 -4*(a*e +b*d)))/2

        else
            % PERIODIC IN Y ONLY

            Ux = KTEMP{3,4};
            Uy = KTEMP{2,3};

			% BLOCKS IN TOP ROW OF LARGE BLOCK FROM A ROW OF ELEMENTS
			D{lvl,1} = KTEMP{3,3}; %first guess
			for nu = 1:1000
				Dnew = KTEMP{3,3} - KTEMP{4,3}*(D{lvl,1}\Uy); 
				if (norm(Dnew -D{lvl,1}) < 1.0e-8) 
					break;
				end
				D{lvl,1} = Dnew;

				if (nu > 999) 
					disp('uh-oh diagonal of top row not converged in ILU');
				end
			end

			Ly{lvl,1} = (KTEMP{4,3}/D{lvl,1})/abs(omega);

			% BLOCKS IN LOWER ROWS OF LARGE BLOCK FROM A ROW OF ELEMENTS

			for N = 1:Nel-1
			    Lx{lvl,N} = (KTEMP{3,2}/D{lvl,N});
			    D{lvl,N+1} = D{lvl,N}; % first guess
			
			    for nu = 1:50000
                    Dnew = KTEMP{3,3}- Lx{lvl,N}*Ux - KTEMP{4,3}*(D{lvl,N+1}\Uy); 
                    if (norm(Dnew -D{lvl,N+1}) < 1.0e-8) 
                        break;
                    end
                    D{lvl,N+1} = Dnew;

                    if (nu > 49999) 
                        disp('uh-oh diagonal of lower rows not converged in ILU');
                    end
                end
                Lx{lvl,N} = Lx{lvl,N}/abs(omega);
                Ly{lvl,N+1} = (KTEMP{4,3}/D{lvl,N+1})/abs(omega);
            end
            
        end
end

% DETERMINE SCALING FOR DIAGONAL BLOCKS
if (~implicit_flag)
    if (omega > 0.0)   
        Nsize = 5;
        if (Fourier2D || oned)
            % For oned Nely = 1 & y-matrices are 0
            % FIND MAXEIG OF RELAXATION SCHEME WITHOUT SWEEPING TERMS
            Kbig = circulant2d(zeros(size(KB{lvl,3,3})),0,0,Nsize,Nsize);
            for m=1:5
                for n=1:5
                    Kbig = Kbig + circulant2d(KB{lvl,m,n},n-3,m-3,Nsize,Nsize);
                end
            end
            Mbig = zeros(size(Kbig));
            for k=hblocks(1):hblocks(2)
                Mbig = Mbig +circulant2d(MB{lvl,3,k},k-3,0,Nsize,Nsize);
            end
        else
            % Use Nelx = 1 in circulant2d
            Kbig = circulant2d(blktimes(zeros(Nsize,Nsize),KB{lvl,3,3}),0,0,1,Nsize);
            for m=1:5
                for n=1:5
                    posy = 3-m;
                    posx = 3-n;
                    Kbig = Kbig + circulant2d(blktimes(diag(ones(Nsize-abs(posx),1),posx),KB{lvl,m,n}),0,posy,1,Nsize);
                end
            end

            Mbig = zeros(size(Kbig));
            for k=hblocks(1):hblocks(2)
                posx = 3-k;
                posy = 0;
                Mbig = Mbig + circulant2d(blktimes(diag(ones(Nsize-abs(posx),1),posx),MB{lvl,3,k}),0,posy,1,Nsize);
            end
        end
        OPTS.disp = 0;
        OPTS.tol = 1.0e-3;
        % biglam = eigs(-Kbig,Mbig,3,'LM',OPTS);
        % biglam = eig(-Kbig,Mbig);
        biglam = eigs(-Mbig\Kbig,3,'LM',OPTS);
        maxeig = max(abs(biglam))
        if (isinf(maxeig))
            biglam = sort(biglam);
            count = size(biglam,1);
            while (isinf(biglam(count)))
                count = count -1;
            end
            biglam = biglam(1:count);
            maxeig = abs(biglam(count));
        end
    else
        maxeig = 1.0;
    end
    
    % SPECIAL CASE!!!
    % REDUCE RELAXATION FACTOR IN DIFFUSION MULTIGRID FOR P=1
    % FOR DIFFUSION AT P=1, BLOCK JACOBI, WITH NO SWEEPING
    if (P == 1  && sys_flag == 0  && (rlx_flag == 2 || rlx_flag == 4) && sweep_flag == 0) 
        omega = 2/3*omega
    end    
    for m=1:5
        for n=1:5
            MB{lvl,m,n} = MB{lvl,m,n}/abs(omega)*maxeig;
            MB1{lvl,m,n} = MB1{lvl,m,n}/abs(omega)*maxeig;
        end
    end
end

% SWEEPING BLOCKS OF RELAXATION MATRIX
if (rlx_flag < 4) 
    switch(sweep_flag)
        case 1
            % GAUSS SEIDEL TYPE
            for m=1:3
                for n=1:3
                    if (m == 3 && n == 3) 
                        continue;
                    end
                 MB{lvl,m,n} = KB{lvl,m,n};
                 MB{lvl,m,n} = KB{lvl,m,n}; 
                end
            end

         case 2
            % SYMMETRIC GAUSS SEIDEL TYPE
            % ASSUME DIAGONAL IS SAME ON BOTH SWEEPS
            MB1{lvl,3,3} = MB{lvl,3,3};

            for m=1:3
                for n=1:3
                    if (m == 3 && n == 3) 
                        continue;
                    end
                 MB{lvl,m,n} = KB{lvl,m,n};
                 MB{lvl,m,n} = KB{lvl,m,n}; 
                end
            end
            
            
            for m=3:5
                for n=3:5
                    if (m == 3 && n == 3) 
                        continue;
                    end
                 MB1{lvl,m,n} = KB{lvl,m,n};
                 MB1{lvl,m,n} = KB{lvl,m,n}; 
                end
            end
    end
end

if (rlx_flag == 4)
    switch(sweep_flag)
        case 1
            % SEQUENTIAL LINES
            % ADD BOTTOM PART
            for m=1:2
                for n=1:5
                 MB{lvl,m,n} = KB{lvl,m,n};
                 MB{lvl,m,n} = KB{lvl,m,n}; 
                end
            end

         case 2
            % SYMMETRIC SEQUENTIAL LINES
            % ADD BOTTOM PART
            for m=1:2
                for n=1:5
                 MB{lvl,m,n} = KB{lvl,m,n};
                 MB{lvl,m,n} = KB{lvl,m,n}; 
                end
            end
           
            % ASSUME LINE IS THE SAME
            for n=1:5
                MB1{lvl,3,n} = MB{lvl,3,n};
            end
            
            % ADD TOP PART
            for m=4:5
                for n=1:5
                 MB1{lvl,m,n} = KB{lvl,m,n};
                 MB1{lvl,m,n} = KB{lvl,m,n}; 
                end
            end
    end
end

if (rlx_flag == 5 && sweep_flag == 1)
    % ALTERNATING DIRECTION LINES WITH SWEEPING
    % ADD BOTTOM PART
    for m=1:2
        for n=1:5
         MB{lvl,m,n} = KB{lvl,m,n};
         MB{lvl,m,n} = KB{lvl,m,n}; 
        end
    end 
    
    % ADD RIGHT PART
    for m=1:5
        for n=1:2
         MB1{lvl,m,n} = KB{lvl,m,n};
         MB1{lvl,m,n} = KB{lvl,m,n}; 
        end
    end 
end
    

if (implicit_flag && (rlx_flag ~= 6 || rlx_flag ~= 7))
    % IMPLICT SCHEME 
    for m=1:5
        for n=1:5
            MB{lvl,m,n} = MB{lvl,m,n}*mu;
            MB1{lvl,m,n} = MB1{lvl,m,n}*mu;
        end
    end
    MB{lvl,3,3} = MB{lvl,3,3} +blktimes(eye(nvar,nvar),ms2d)*dx*dy/dt;
    MB1{lvl,3,3} = MB1{lvl,3,3} +blktimes(eye(nvar,nvar),ms2d)*dx*dy/dt;
end


% THIS IS TO DO AGGLOMERATION MULTIGRID
if (P==0 && lvl==1 && maxlvl > 1)
    disp('Agglomerating');
    % Elements are order counterclockwise from bottom left
    if (oned)
        % THIS DOES 2 1D LINES AND ONLY RESTRICTS IN X DIRECTION
        restrict2d = [1 1 0 0; 0 0 1 1];
    else
        restrict2d= [1 1 1 1];
    end

    %   REARRANGEMENT OF THE STIFFNESS MATRIX
    KB{lvl,3,3} = blktimes(KB{lvl,3,3},eye(4,4));
    KB{lvl,3,3} = KB{lvl,3,3} + blktimes(KB{lvl,4,3},[0,0,0,1; 0,0,1,0; 0,0,0,0; 0,0,0,0]);
    KB{lvl,3,3} = KB{lvl,3,3} + blktimes(KB{lvl,2,3},[0,0,0,0; 0,0,0,0; 0,1,0,0; 1,0,0,0]);
    KB{lvl,3,3} = KB{lvl,3,3} + blktimes(KB{lvl,3,2},[0,0,0,0; 1,0,0,0; 0,0,0,1; 0,0,0,0]);
    KB{lvl,3,3} = KB{lvl,3,3} + blktimes(KB{lvl,3,4},[0,1,0,0; 0,0,0,0; 0,0,0,0; 0,0,1,0]);
    
    KB{lvl,3,4} = blktimes(KB{lvl,3,4},[0,0,0,0; 1,0,0,0; 0,0,0,1; 0,0,0,0]);
    KB{lvl,4,3} = blktimes(KB{lvl,4,3},[0,0,0,0; 0,0,0,0; 0,1,0,0; 1,0,0,0]);
    
    KB{lvl,3,2} = blktimes(KB{lvl,3,2},[0,1,0,0; 0,0,0,0; 0,0,0,0; 0,0,1,0]);
    KB{lvl,2,3} = blktimes(KB{lvl,2,3},[0,0,0,1; 0,0,1,0; 0,0,0,0; 0,0,0,0]);
    KB{lvl,3,1}= zeros(size(KB{lvl,3,3})); KB{lvl,3,5}= zeros(size(KB{lvl,3,3})); KB{lvl,5,3}= zeros(size(KB{lvl,3,3})); KB{lvl,1,3}= zeros(size(KB{lvl,3,3}));
    KB{lvl,1,1}= zeros(size(KB{lvl,3,3})); KB{lvl,2,2}= zeros(size(KB{lvl,3,3})); KB{lvl,1,2}= zeros(size(KB{lvl,3,3})); KB{lvl,2,1}= zeros(size(KB{lvl,3,3}));
    KB{lvl,4,1}= zeros(size(KB{lvl,3,3})); KB{lvl,4,2}= zeros(size(KB{lvl,3,3})); KB{lvl,5,1}= zeros(size(KB{lvl,3,3})); KB{lvl,5,2}= zeros(size(KB{lvl,3,3}));
    KB{lvl,1,4}= zeros(size(KB{lvl,3,3})); KB{lvl,1,5}= zeros(size(KB{lvl,3,3})); KB{lvl,2,4}= zeros(size(KB{lvl,3,3})); KB{lvl,2,5}= zeros(size(KB{lvl,3,3}));
    KB{lvl,4,4}= zeros(size(KB{lvl,3,3})); KB{lvl,4,5}= zeros(size(KB{lvl,3,3})); KB{lvl,5,4}= zeros(size(KB{lvl,3,3})); KB{lvl,5,5}= zeros(size(KB{lvl,3,3}));
    
    %   REARRANGEMENT OF THE RELAXATION MATRIX
    MB{lvl,3,3} = blktimes(MB{lvl,3,3},eye(4,4));
    MB{lvl,3,3} = MB{lvl,3,3} + blktimes(MB{lvl,4,3},[0,0,0,1; 0,0,1,0; 0,0,0,0;0,0,0,0]);
    MB{lvl,3,3} = MB{lvl,3,3} + blktimes(MB{lvl,2,3},[0,0,0,0; 0,0,0,0; 0,1,0,0; 1,0,0,0]);
    MB{lvl,3,3} = MB{lvl,3,3} + blktimes(MB{lvl,3,2},[0,0,0,0; 1,0,0,0; 0,0,0,1; 0,0,0,0]);
    MB{lvl,3,3} = MB{lvl,3,3} + blktimes(MB{lvl,3,4},[0,1,0,0; 0,0,0,0; 0,0,0,0; 0,0,1,0]);
    
    MB{lvl,3,4} = blktimes(MB{lvl,3,4},[0,0,0,0; 1,0,0,0; 0,0,0,1; 0,0,0,0]);
    MB{lvl,4,3} = blktimes(MB{lvl,4,3},[0,0,0,0; 0,0,0,0; 0,1,0,0; 1,0,0,0]);
    
    MB{lvl,3,2} = blktimes(MB{lvl,3,2},[0,1,0,0; 0,0,0,0; 0,0,0,0; 0,0,1,0]);
    MB{lvl,2,3} = blktimes(MB{lvl,2,3},[0,0,0,1; 0,0,1,0; 0,0,0,0; 0,0,0,0]);
    MB{lvl,3,1}= zeros(size(MB{lvl,3,3})); MB{lvl,3,5}= zeros(size(MB{lvl,3,3})); MB{lvl,5,3}= zeros(size(MB{lvl,3,3})); MB{lvl,1,3}= zeros(size(MB{lvl,3,3}));
    MB{lvl,1,1}= zeros(size(KB{lvl,3,3})); MB{lvl,2,2}= zeros(size(KB{lvl,3,3})); MB{lvl,1,2}= zeros(size(KB{lvl,3,3})); MB{lvl,2,1}= zeros(size(KB{lvl,3,3}));
    MB{lvl,4,1}= zeros(size(KB{lvl,3,3})); MB{lvl,4,2}= zeros(size(KB{lvl,3,3})); MB{lvl,5,1}= zeros(size(KB{lvl,3,3})); MB{lvl,5,2}= zeros(size(KB{lvl,3,3}));
    MB{lvl,1,4}= zeros(size(KB{lvl,3,3})); MB{lvl,1,5}= zeros(size(KB{lvl,3,3})); MB{lvl,2,4}= zeros(size(KB{lvl,3,3})); MB{lvl,2,5}= zeros(size(KB{lvl,3,3}));
    MB{lvl,4,4}= zeros(size(KB{lvl,3,3})); MB{lvl,4,5}= zeros(size(KB{lvl,3,3})); MB{lvl,5,4}= zeros(size(KB{lvl,3,3})); MB{lvl,5,5}= zeros(size(KB{lvl,3,3}));
    
    
    % REARRANGEMENT OF THE RELAXATION MATRIX FOR SYMMETRIC SWEEPS
    MB1{lvl,3,3} = blktimes(MB1{lvl,3,3},eye(4,4));
    MB1{lvl,3,3} = MB1{lvl,3,3} + blktimes(MB1{lvl,4,3},[0,0,0,1; 0,0,1,0; 0,0,0,0;0,0,0,0]);
    MB1{lvl,3,3} = MB1{lvl,3,3} + blktimes(MB1{lvl,2,3},[0,0,0,0; 0,0,0,0; 0,1,0,0; 1,0,0,0]);
    MB1{lvl,3,3} = MB1{lvl,3,3} + blktimes(MB1{lvl,3,2},[0,0,0,0; 1,0,0,0; 0,0,0,1; 0,0,0,0]);
    MB1{lvl,3,3} = MB1{lvl,3,3} + blktimes(MB1{lvl,3,4},[0,1,0,0; 0,0,0,0; 0,0,0,0; 0,0,1,0]);
    
    MB1{lvl,3,4} = blktimes(MB1{lvl,3,4},[0,0,0,0; 1,0,0,0; 0,0,0,1; 0,0,0,0]);
    MB1{lvl,4,3} = blktimes(MB1{lvl,4,3},[0,0,0,0; 0,0,0,0; 0,1,0,0; 1,0,0,0]);
    
    MB1{lvl,3,2} = blktimes(MB1{lvl,3,2},[0,1,0,0; 0,0,0,0; 0,0,0,0; 0,0,1,0]);
    MB1{lvl,2,3} = blktimes(MB1{lvl,2,3},[0,0,0,1; 0,0,1,0; 0,0,0,0; 0,0,0,0]);
    MB1{lvl,3,1}= zeros(size(MB1{lvl,3,3})); MB1{lvl,3,5}= zeros(size(MB1{lvl,3,3})); MB1{lvl,5,3}= zeros(size(MB1{lvl,3,3})); MB1{lvl,1,3}= zeros(size(MB1{lvl,3,3}));
    
    MB1{lvl,1,1}= zeros(size(KB{lvl,3,3})); MB1{lvl,2,2}= zeros(size(KB{lvl,3,3})); MB1{lvl,1,2}= zeros(size(KB{lvl,3,3})); MB1{lvl,2,1}= zeros(size(KB{lvl,3,3}));
    MB1{lvl,4,1}= zeros(size(KB{lvl,3,3})); MB1{lvl,4,2}= zeros(size(KB{lvl,3,3})); MB1{lvl,5,1}= zeros(size(KB{lvl,3,3})); MB1{lvl,5,2}= zeros(size(KB{lvl,3,3}));
    MB1{lvl,1,4}= zeros(size(KB{lvl,3,3})); MB1{lvl,1,5}= zeros(size(KB{lvl,3,3})); MB1{lvl,2,4}= zeros(size(KB{lvl,3,3})); MB1{lvl,2,5}= zeros(size(KB{lvl,3,3}));
    MB1{lvl,4,4}= zeros(size(KB{lvl,3,3})); MB1{lvl,4,5}= zeros(size(KB{lvl,3,3})); MB1{lvl,5,4}= zeros(size(KB{lvl,3,3})); MB1{lvl,5,5}= zeros(size(KB{lvl,3,3}));
end    

% MAKE RESTRICTION MATRIX (USUALLY ONLY DIAGONAL IS NON-ZERO EXCEPT IN RESTRICTION TO C0)
for m=1:3
    for n=1:3
        RB{lvl,n,m} = zeros(nvar*size(restrict2d));
    end
end
RB{lvl,2,2} = blktimes(eye(nvar,nvar),restrict2d);

if (cgswitch && lvl == maxlvl-1) 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FOLLOWING LINES ARE FOR COARSENING TO CONTINUOUS SPACE AFTER P = 1 %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % harold_tri_test;
    
    if (basis_flag ~= 4)
        % ASSUMES NODAL BASIS WITH ORDERING 
        % 3-4
        % | |
        % 1-2
        disp('Must use nodal basis for restriction to continuous space');
    end
    if (~supg_restrict)
        if (oned)
            RB{lvl,2,2} = blktimes(eye(nvar,nvar),[1,0]);
            RB{lvl,2,1} = blktimes(eye(nvar,nvar),[0,1]);
        else
            RB{lvl,2,2} = blktimes(eye(nvar,nvar),[1,0,0,0]);
            RB{lvl,2,1} = blktimes(eye(nvar,nvar),[0,1,0,0]);
            RB{lvl,1,2} = blktimes(eye(nvar,nvar),[0,0,1,0]);
            RB{lvl,1,1} = blktimes(eye(nvar,nvar),[0,0,0,1]);
        end
    else
        if (oned) 
            % This constructs V +ax*tau*dV/dxi on coarse mesh given DG residuals
            % dv/dx of test function negative constant = 1/2
            RB{lvl,2,2} = blktimes(eye(nvar,nvar),[1,0]) +blktimes(ax*tau,[-0.5,-0.5]);
            % dv/dx of test function positive constant = 1/2
            RB{lvl,2,2} = blktimes(eye(nvar,nvar),[1,0]) +blktimes(ax*tau,[0.5,0.5]);
        else
            % This constructs V +ax*tau*dV/dxi +ay*tau*dV/deta on coarse mesh
            % dv/dx dv/dy negative
            RB{lvl,2,2} = blktimes(eye(nvar,nvar),[1,0,0,0]) +blktimes(ax*tau,[-0.5,-0.5,0,0]) +blktimes(ay*tau,[-0.5,0,-0.5,0]);
            % dv/dx positive dv/dy negative
            RB{lvl,2,1} = blktimes(eye(nvar,nvar),[0,1,0,0]) +blktimes(ax*tau,[0.5,0.5,0,0]) +blktimes(ay*tau,[-0.5,0,-0.5,0]);
            % dv/dx negative dv/dy positive
            RB{lvl,1,2} = blktimes(eye(nvar,nvar),[0,0,1,0]) +blktimes(ax*tau,[0,0,-0.5,-0.5]) +blktimes(ay*tau,[0.5,0,0.5,0]);
            % dv/dx positive dv/dy positive
            RB{lvl,1,1} = blktimes(eye(nvar,nvar),[0,0,0,1]) +blktimes(ax*tau,[0.0,0.0,0.5,0.5]) +blktimes(ay*tau,[0.0,0.5,0.0,0.5]);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TESTING TO CHECK 2-LEVEL RESULTS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_flag = 0;
if (lvl == 1 && test_flag == 1)

    if (Fourier2D || oned)
        % For oned Nely = 1 & y-matrices are 0
        Kbig = circulant2d(zeros(size(KB{lvl,3,3})),0,0,Nelx,Nely);
        Mbig = zeros(size(Kbig));
        Mbig1 = zeros(size(Kbig));
        Rbig = circulant2d(zeros(size(RB{lvl,2,2})),0,0,Nelx,Nely);
        for m=1:5
            for n=1:5
                Kbig = Kbig + circulant2d(KB{lvl,m,n},n-3,m-3,Nelx,Nely);
                Mbig = Mbig + circulant2d(MB{lvl,m,n},n-3,m-3,Nelx,Nely);
                Mbig1 = Mbig1 + circulant2d(MB1{lvl,m,n},n-3,m-3,Nelx,Nely);
            end
        end

        for m=1:3
            for n=1:3
                Rbig = Rbig + circulant2d(RB{lvl,m,n},n-2,m-2,Nelx,Nely);
            end
        end
    else
        % Use Nelx = 1 in circulant2d
        Kbig = circulant2d(blktimes(zeros(Nelx,Nelx),KB{lvl,3,3},0,0,1,Nely));
        Mbig = zeros(size(Kbig));
        Mbig1 = zeros(size(Kbig));
        Rbig = circulant2d(blktimes(zeros(Nelx,Nelx),RB{lvl,2,2}),0,0,1,Nely);
        for m=1:5
            for n=1:5
                posy = 3-m;
                posx = 3-n;
                Kbig = Kbig + circulant2d(blktimes(diag(ones(Nelx-abs(posx),1),posx),KB{lvl,m,n}),0,posy,1,Nely);
                Mbig = Mbig + circulant2d(blktimes(diag(ones(Nelx-abs(posx),1),posx),MB{lvl,m,n}),0,posy,1,Nely);
                Mbig1 = Mbig1 + circulant2d(blktimes(diag(ones(Nelx-abs(posx),1),posx),MB1{lvl,m,n}),0,posy,1,Nely);
            end
        end

        for m=1:2
            for n=1:2
                posy = 2-m;
                posx = 2-n;
                Rbig = Rbig + circulant2d(blktimes(diag(ones(Nelx-abs(posx),1),posx),RB{lvl,m,n}),0,posy,1,Nely);
            end
        end
    end
        

    % TESTING RELAXATION SCHEME
    fine = (eye(size(Kbig)) -inv(Mbig)*Kbig);
    if (sweep_flag == 2)
        % SECOND STEP OF SYMMETRIC SWEEP
        fine = (eye(size(Kbig)) -inv(Mbig1)*Kbig)*fine;
    end
    fineeig = eig(fine);
    figure
    plot(real(fineeig),imag(fineeig),'x');
    title('amplification factor (direct)');
    axis equal
    ayactordirect =max(abs(fineeig(2:size(fineeig,1))))
    
    % TESTING TWO-LEVEL MULTIGRID SCHEME
    Kbigc = Rbig*Kbig*(Rbig');
    wholething = (eye(size(Kbig)) -Rbig'*lscov(Kbigc,Rbig*Kbig))*fine;
%   wholething = (eye(size(Kbig)) -Rbig'*(Kbigc\(Rbig*Kbig)))*fine;
    eigswt =eig(wholething);
    figure
    plot(real(eigswt),imag(eigswt),'x')
    title('multigrid amplification factor (direct)');
    axis equal
    dfactordirect =max(abs(eigswt(2:size(eigswt,1))))
 end



