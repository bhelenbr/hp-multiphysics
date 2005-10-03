% RELOAD FINE GRID MATRICES
lvl = 1;

maxeig = vmaxeig(lvl);
bsz = size(KB{lvl,3,3},1);


NDIVX = 18;
NDIVY = NDIVX;

global kdx kdy;

% KDX MISSES 0 
% KDY HITS 0
for kcntx=1:2*NDIVX
    kdx = (kcntx-NDIVX-1/2)*pi/NDIVX;
    
	for kcnty=1:2*NDIVY-1
        kdy = (kcnty-NDIVY)*pi/NDIVY;
        
%         kdx = 0.0;
%         kdy = pi;
        
		% STIFFNESS
        k1 = zeros(bsz,bsz);
        for m=1:5
            for n=1:5
                k1 = k1 +KB{lvl,m,n}*exp(i*((m-3)*kdy +(n-3)*kdx));
            end
        end
        
        % PRECONDITIONER
        if (sweep_flag ~= 2 && sweep_flag ~= 5)
            m1 = zeros(bsz,bsz);
            for m=1:5
                for n=1:5
                    m1 = m1 +MB{lvl,m,n}*exp(i*((m-3)*kdy +(n-3)*kdx));
                end
            end
            m1 = inv(m1);
        elseif (sweep_flag == 2)
            mbl = zeros(bsz,bsz);
            for m=1:3
                for n=1:3
                    mbl = mbl +MB{lvl,m,n}*exp(i*((m-3)*kdy +(n-3)*kdx));
                end
            end
            mtr = zeros(bsz,bsz);
            for m=3:5
                for n=3:5
                    mtr = mtr +MB{lvl,m,n}*exp(i*((m-3)*kdy +(n-3)*kdx));
                end
            end   
            m1 = inv(mtr) +inv(mbl) -inv(mtr)*k1*inv(mbl);
        else
            mb = zeros(bsz,bsz);
            for m=1:3
                for n=1:5
                    mb = mb +MB{lvl,m,n}*exp(i*((m-3)*kdy +(n-3)*kdx));
                end
            end
            mt = zeros(bsz,bsz);
            for m=3:5
                for n=1:5
                    mt = mt +MB{lvl,m,n}*exp(i*((m-3)*kdy +(n-3)*kdx));
                end
            end   
            m1 = inv(mt) +inv(mb) -inv(mt)*k1*inv(mb);
        end
            
            
        % FINE GRID EIGENVALUES
        a = -m1*k1;
        ea = eig(a);
        ea = sort(ea);
        
        % FINE GRID AMPLIFICATION FACTOR
        if (sys_flag == 1 || rlx_flag == 2)
            rka = eye(size(a)) +a;
        else
            rka = rk3_5(eye(size(a)),zeros(size(a)),-2.0*a,zeros(size(a)));
        end
        erka = eig(rka);
        [tmp,ind] = sort(abs(abs(erka)-1.0));
        erka = erka(ind);
        
		lam(kcntx,kcnty,:) = [ea.'];
        amp(kcntx,kcnty,:) = [erka.'];
            
        %%%%%%%%%%%%%%%%%%%%
        % MULTIGRID ANALYSIS
        %%%%%%%%%%%%%%%%%%%%
        if (maxlvl > 1)
            mgamp1 = mgcycle2d(vw,1,0);
         
            mgamp = eig(mgamp1);
            mgamp = sort(mgamp);
            
			mgrid(kcntx,kcnty,:) = [mgamp.'];
        end
        
        xpos(kcntx,1:2*NDIVY-1) = kdx;
        ypos(1:2*NDIVX,kcnty) = kdy;
        
%         %%%% TO PLOT EIGENFUNCTIONS %%%%
%         [vmg,emg] = eig(a);
%         emg = diag(emg);
%         [emg,smg] = sort(emg);
%         vmg = vmg(:,smg);
%         
%         for vct=1:bsz;
%             figure;
%             elmnt = 0;
%             hold on;
%             mode = basis2dsave*vmg(:,vct);
%             syms lx ly
%             if (kdx == 0) 
%                 NX = 1;
%             else
%                 NX = 2*pi/kdx;
%             end
%             for ex=0:NX-1
%                 for ey=0:2*pi/kdy-1
%                     ezcontourf(real(subs(mode*exp(i*(ex*kdx +ey*kdy)),{x,y},{lx-ex*2,ly-ey*2})),[-1+ex*2,1+ex*2,-1+ey*2,1+ey*2]);
%                 end
%             end
%             axis auto;
%             axis equal;
%             title(num2str(emg(vct)));
%         end
%         return;
    end
end

%xcut = NDIVY;
%disp('ZONE');
%bigcat(:,1) = xpos(:,1);
%for ic = 1:bsz
 %   bigcat(:,ic+1) = 1.0-abs(amp(:,1,ic));
 %end
%disp(num2str(bigcat));
%return;

figure;
for ic=1:bsz
plot(real(lam(:,:,ic)),imag(lam(:,:,ic)),'x');
hold on;
end
axis equal;

sqrtbsz = sqrt(bsz);

% figure;
% hold off;
% for ic=1:bsz
% subplot(sqrtbsz,sqrtbsz,ic);
% contourf(xpos,ypos,abs(lam(:,:,ic)));
% colorbar;
% end
% title('\lambda');
% % 
% figure;
% for ic=1:bsz
% 	subplot(sqrtbsz,sqrtbsz,ic);
% 	contourf(xpos,ypos,abs(amp(:,:,ic)));
% 	colorbar;
% end
% title('Amplification Factor');
rfactor = max(max(max(abs(amp(:,:,:)))))

if (maxlvl > 1)
	figure;
% %    system('rm mgrid.dat');
% 	for ic=1:bsz
% 	subplot(sqrtbsz,sqrtbsz,ic);
% 	contourf(xpos,ypos,abs(mgrid(:,:,ic)));
% 	colorbar;
% %    matrix2tecplot(xpos,ypos,abs(mgrid(:,:,bsz+1-ic)),'mgrid.dat');
% 	end
% 	title('Multigrid Damping Factor');
	[dfactor,IP] = max(abs(mgrid(:,:,bsz)));
    [dfactor,JP] = max(dfactor);
    IP = IP(JP);
    dfactor
    kdxmax = xpos(IP,JP)
    kdymax = ypos(IP,JP)
end


% CHECK 1-D SLICE AT K DY = 0
xcut = NDIVY;

figure;
hold off;
subplot(2,2,1);
for ic=1:bsz
plot(xpos(:,xcut),abs(lam(:,xcut,ic)),'+');
hold on;
end
ylabel('|\lambda|','FontSize',14);
xlabel('k dx','FontSize',14);

hold off;
subplot(2,2,2);
for ic=1:bsz
    plot(real(lam(:,xcut,ic)),imag(lam(:,xcut,ic)),'+');
    hold on;
end
axis equal;
xlabel('Real(\lambda)','FontSize',14);
ylabel('Imag(\lambda)','FontSize',14);

hold off;
subplot(2,2,3);
for ic=1:bsz
plot(xpos(:,xcut),abs(amp(:,xcut,ic)),'+');
hold on;
end
ylabel('Fine Grid Amplification Factor','FontSize',14);
xlabel('k dx','FontSize',14);

if (maxlvl > 1)
	hold off;
	subplot(2,2,4);
	for ic=1:bsz
	plot(xpos(:,xcut),abs(mgrid(:,xcut,ic)),'+');
	hold on;
	end
	ylabel('Multigrid Amplification Factor','FontSize',14);
	xlabel('k dx','FontSize',14);
    % dfactor1d = max(abs(mgrid(:,xcut,:)))
end