% RELOAD FINE GRID MATRICES
lvl = 1;
maxeig = vmaxeig{lvl};
restrict = vrestrict{lvl};
bsz = size(K,1);

global kdx;

NDIV = 49; 

plotfuncs = 0
if (plotfuncs) 
    NDIV = 4;
    kdxarray = 2*pi./[16 8 4 2]
end

for kcnt=1:NDIV
	kdx = kcnt*pi/(NDIV+1);
    
    if (plotfuncs)
        kdx = kdxarray(kcnt);
    end
          
	% STIFFNESS
    k1 = zeros(size(KB{lvl,3}));
    for m = 1:5
	    k1 = k1 +KB{lvl,m}*exp(-i*(m-3)*kdx);
    end
    
	% RELAXATION
    if (sweep_flag ~= 2) 
		m1 = zeros(size(KB{lvl,3}));
		for m = 1:5
            m1 = m1 + MB{lvl,m}*exp(-i*(m-3)*kdx);
		end
        a = -inv(m1)*k1;
    else
        ml = zeros(size(KB{lvl,3}));
		for m = 1:3
            ml = ml + MB{lvl,m}*exp(-i*(m-3)*kdx);
		end   
        mr = zeros(size(KB{lvl,3}));
		for m = 3:5
            mr = mr + MB{lvl,m}*exp(-i*(m-3)*kdx);
		end   
        m1 = inv(mr) +inv(ml) -inv(mr)*k1*inv(ml);
        a = -m1*k1;
    end

    
    % FINE GRID EIGENVALUES
    ea = eig(a);
    [tmp,ind] = sort(abs(ea));
    ea = ea(ind);
    
    % FINE GRID AMPLIFICATION FACTOR
    if (sys_flag == 1 || rlx_flag == 2)
        rka = eye(size(a)) +a;
    else
        % CFL OF 2:
        rka = rk3_5(eye(size(a)),zeros(size(a)),-2.0*a,zeros(size(a)));
    end
    erka = eig(rka);
    [tmp,ind] = sort(abs(abs(erka)-1.0));
    erka = erka(ind);
    
	lam(kcnt,:) = [kdx, ea.'];
    amp(kcnt,:) = [kdx, erka.'];
	
	%%%%%%%%%%%%%%%%%%%%
    % MULTIGRID ANALYSIS
    %%%%%%%%%%%%%%%%%%%%
    mgamp1 = mgcycle(vw,1,0);
 
    mgamp = eig(mgamp1);
    mgamp = sort(mgamp);
    
	mgrid(kcnt,:) = [kdx, mgamp.'];
    

%   TO PLOT EIGENFUNCTIONS
    if (plotfuncs) 
        [vmg,emg] = eig(rka);
        emg = diag(emg);
        [emg,smg] = sort(emg);
        vmg = vmg(:,smg);
        
        [vmga,emga] = eig(mgamp1);
        emga = diag(emga);
        [emga,smga] = sort(emga);
        vmga = vmga(:,smga);
        
        
        figure;

		for vct = 1:size(vmg,1);
            subplot(size(vmg,1),3,3*vct-2);
%            subplot(size(vmg,1),1,vct);
            hold on;
			mode = basissave*vmg(:,vct);
			syms lx
			if (kdx == 0) 
                NX = 1;
			else
                NX = 2*pi/kdx;
			end
			for ex=0:NX-1
                ezplot(real(subs(mode*exp(-i*ex*kdx),x,(2*(lx-1/2)-ex*2))),1/2 +[-1+ex*2,1+ex*2]/2);
                xlabel('x/\Deltax');
			end
            axis auto;            
            title(['damping factor:' num2str(abs(emg(vct))) ' theta/pi:' num2str(kdx/pi)] );
            
            subplot(size(vmg,1),3,3*vct-1);
            hold on;
			mode = basissave*vmga(:,vct);
			syms lx
			if (kdx == 0) 
                NX = 1;
			else
                NX = 2*pi/kdx;
			end
			for ex=0:NX-1
                ezplot(real(subs(mode*exp(-i*ex*kdx),x,lx-ex*2)),[-1+ex*2,1+ex*2]);
			end
            axis auto;
            title(['multigrid damping factor:' num2str(abs(emga(vct))) ' theta/pi:' num2str(kdx/pi)] );
            
            subplot(size(vmg,1),3,3*vct);
            hold on;
            
                    % RESTRICTION OPERATOR
            r1 = zeros(size(RB{lvl,2}));
            for m = 1:3
                r1 = r1 + RB{lvl,m}*exp(-i*(m-2)*kdx);
            end
            dampedmode = vmg(:,vct) -r1'*mgcycle(vw,2,r1*k1*vmg(:,vct));

			mode = basissave*dampedmode;
			syms lx
			if (kdx == 0) 
                NX = 1;
			else
                NX = 2*pi/kdx;
			end
			for ex=0:NX-1
                ezplot(real(subs(mode*exp(-i*ex*kdx),x,lx-ex*2)),[-1+ex*2,1+ex*2]);
			end
            axis auto;
            title(['mgrid on fine function: theta/pi:' num2str(kdx/pi)] );
		end
    end
end

bsz = size(lam,2)-1;
% COMPARISON TO ANALYTIC FOR DIFFUSION
if (sys_flag == 1 & rlx_flag+sweep_flag == 0)
% 	figure;
% 	hold off;
% 	for ic=1:bsz
% 	plot(floor(ic/2)*2*pi +(-1)^(ic-1)*lam(:,1),-4.*real(lam(:,ic+1))./(floor(ic/2)*2*pi +(-1)^(ic-1)*lam(:,1)).^2,'+');
% 	hold on;
% 	end
% 	xlabel('k dx','FontSize',14);
% 	ylabel('- \lambda / k^2','FontSize',14);
% 	%outstring = ['/Users/helenbrk/Codes/analysis/multigrid/DG viscous plots/' fstring '.eps']
% 	%print('-deps',outstring)
	
	% COMPARISON TO ANALYTIC FOR DIFFUSION
    %figure;
    if (scheme == 7) 
        scheme = 4
    end
    strings1 = {'bo-','gx-','r+-','m*-','ks-','bd-'};
    strings2 = {'Bassi-Rebay','Brezzi et al.','local DG \beta=0','local DG \beta=1/2, \eta=0','interior penalty','Bassi et al'};
	for ic=1:bsz
    loglog(floor(ic/2)*2*pi +(-1)^(ic-1)*lam(:,1),abs(-4.*real(lam(:,ic+1))-(floor(ic/2)*2*pi +(-1)^(ic-1)*lam(:,1)).^2),char(strings1(scheme)));

    hold on;
	end
    plot(1.2*10^-1,100.0/10^(scheme/1.1),char(strings1(scheme)));
    text(1.4*10^-1,100.0/10^(scheme/1.1),char(strings2(scheme)),'FontSize',14);
    axis([0.1 floor(bsz/2)*2*pi+pi 10^(-14) 10^4])
	xlabel('\theta','FontSize',14);
	ylabel('error in eigenvalue','FontSize',14);
	%outstring = ['/Users/helenbrk/Codes/analysis/multigrid/DG viscous plots/' fstring '.eps']
	%print('-deps',outstring)
elseif (sys_flag == 0 & rlx_flag+sweep_flag == 0)
    % COMPARISON TO ANALYTIC FOR CONVECTION
	figure;
	hold off;
	for ic=1:bsz
	plot(floor(ic/2)*2*pi +(-1)^(ic-1)*lam(:,1),abs(2.*imag(lam(:,ic+1)))./(floor(ic/2)*2*pi +(-1)^(ic-1)*lam(:,1)),'+');
	hold on;
	end
	xlabel('k dx','FontSize',14);
	ylabel('|\lambda| / k','FontSize',14);
	%outstring = ['/Users/helenbrk/Codes/analysis/multigrid/DG viscous plots/' fstring '.eps']
	%print('-deps',outstring)
	
	% COMPARISON TO ANALYTIC FOR CONVECTION
    figure;
	for ic=1:bsz
	loglog(floor(ic/2)*2*pi +(-1)^(ic-1)*lam(:,1),abs(abs(-2.*imag(lam(:,ic+1)))./(floor(ic/2)*2*pi +(-1)^(ic-1)*lam(:,1))-1),'m');
	hold on;
	end
	xlabel('k dx','FontSize',14);
	ylabel('Error','FontSize',14);
	%outstring = ['/Users/helenbrk/Codes/analysis/multigrid/DG viscous plots/' fstring '.eps']
	%print('-deps',outstring)
end
    
% PLOT MAGNITUDE OF EIGENVALUES
figure;
hold off;
%subplot(2,2,1)
for ic=1:bsz
plot(floor(ic/2)*2*pi +(-1)^(ic-1)*lam(:,1),abs(lam(:,ic+1)),'+');
hold on;
end
ylabel('|\lambda|','FontSize',14);
xlabel('k dx','FontSize',14);

figure;
hold off;
%subplot(2,2,2);
for ic=1:bsz
plot(real(lam(:,ic+1)),abs(imag(lam(:,ic+1))),'+');
hold on;
end
axis equal;
xlabel('Real(\lambda)','FontSize',14);
ylabel('Imag(\lambda)','FontSize',14);

% figure;
% hold off;
% %subplot(2,2,2);
% for ic=1:bsz
% plot(real(amp(:,ic+1)),abs(imag(amp(:,ic+1))),'+');
% hold on;
% end
% axis equal;
% xlabel('Real(amp)','FontSize',14);
% ylabel('Imag(amp)','FontSize',14);
% 
% 
% figure;
% hold off;
% %subplot(2,2,2);
% for ic=1:bsz
% plot(real(mgrid(:,ic+1)),abs(imag(mgrid(:,ic+1))),'+');
% hold on;
% end
% axis equal;
% xlabel('Real(mgrid)','FontSize',14);
% ylabel('Imag(mgrid)','FontSize',14);

figure;
hold off;
rfactor = 0;
%subplot(2,2,3);
for ic=1:bsz
plot(floor(ic/2)*2*pi +(-1)^(ic-1)*amp(:,1),abs(amp(:,ic+1)),'+');  % TEMPORARY
rfactor = max(max(abs(amp(:,ic+1))),rfactor);
hold on;
end
rfactor
ylabel('damping factor','FontSize',14);
xlabel('\theta','FontSize',14);

figure;
hold off;
%subplot(2,2,4)
dfactor=0;
for ic=1:size(mgrid,2)-1
plot(floor(ic/2)*2*pi +(-1)^(ic-1)*mgrid(:,1),abs(mgrid(:,ic+1)),'gx');
dfactor = max(max(abs(mgrid(:,ic+1))),dfactor);
hold on;
end
dfactor
axis([0 floor(bsz/2)*2*pi+pi 0 max(1.0,dfactor)])

ylabel('multigrid damping factor','FontSize',14);
xlabel('\theta','FontSize',14);

% title(tstring,'FontSize',18);
% xlabel('k dx','FontSize',14);
% ylabel('Multigrid Damping Factor','FontSize',14);
% outstring = ['/Users/helenbrk/Codes/analysis/multigrid/DG viscous plots/' fstring '.mg.eps']
% print('-deps',outstring)
