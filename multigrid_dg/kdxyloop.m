global kdx kdy;
kdxarray = -pi:2*pi/Nelx:pi-2*pi/Nelx;
kdyarray = -pi:2*pi/Nely:pi-2*pi/Nely;

% TO AVOID kdx=0
%kdxarray = -pi+pi/Nelx:2*pi/Nelx:pi

% CHEATS FOR SPEED (4 elements in y-direction)
%kdyarray = -pi:pi/4:pi-pi/4
mgridonly = 0;
fastway = 0;

% SET TO 1 TO PLOT RELAXATION EIGENFUNCTIONS
% SET TO 2 TO PLOT RELAXATION & MULTIGRID EIGENFUNCTIONS
% SET TO 3 TO PLOT RELAXATION & MULTIGRID & MULTIGRID APPLIED TO RELAXATION EIGENFUNCTION
plotfuncs=0;
if (plotfuncs)
    % PICK WAVENUMBER TO PLOT HERE
    kdxarray = 2*pi/Nelx;
    kdyarray = -pi/4; %-pi/4;
    kdxarray = 0;
    howmany = 1;  % 1 or bsz
end

% SET TO 0 FOR NO PLOTTING
% SET TO 1 FOR COMPLEX EIGENVALUE PLOTS
% SET TO 2 FOR ADDITIONAL MORE DETAILED PLOTS
ploteigs=0;
if (oned)
    kdyarray = 0.0;
elseif (~Fourier2D)
    kdxarray = 0.0;
end

% Parameters for faster eigenvalue solves when no plotting is necessary
if (fastway)
    OPTSmg.disp = 0;
    OPTSmg.tol = 1.0e-6;

    if (oned || Fourier2D)
        vecsize = size(KB{1,3,3},1);
    else
        vecsize = Nelx*size(KB{1,3,3},1);
    end
    OPTSmg.maxit = max(vecsize,300);
    OPTSmg.v0 = rand(vecsize,1);

    OPTSrlx = OPTSmg;
end

%FOURIER ANALYSIS 
for kcntx=1:size(kdxarray,2)
    kdx = kdxarray(kcntx);
    
	for kcnty=1:size(kdyarray,2)
        kdy = kdyarray(kcnty);

        xpos(kcntx,kcnty) = kdx;
        ypos(kcntx,kcnty) = kdy;
        
        
        if (abs(kdx) < 1.0e-10 && abs(kdy) < 1.0e-10 && (oned || Fourier2D))
            % AVOID 0 SINGULAR POINT?
            if (ploteigs) 
                % SHIFT A LITTLE SO DON'T PUT A BLANK 0 IN THE PLOTS
                kdx = pi/Nelx;
                xpos(kcntx,kcnty) = kdx;
            else
                % JUST SKIP IT
                amp(kcntx,kcnty,:) = 0.0;
                mgrid(kcntx,kcnty,:) = 0.0;
                continue
            end
        end

        [mgamp1,rka] = mgcycle2d(vw,1,0);
        
        % FINE GRID EIGENVALUES
        if (ploteigs || ~fastway)
            erka = eig(rka);
            [tmp,ind] = sort(abs(abs(erka)-1.0));  % SORTED BASED ON EIGENVALUE NOT DAMPING FACTOR
            erka = erka(ind);
            amp(kcntx,kcnty,:) = [erka.'];
  
            if (maxlvl > 1)
                %%%%%%%%%%%%%%%%%%%%%%%
                % MULTIGRID EIGENVALUES
                %%%%%%%%%%%%%%%%%%%%%%%
                mgamp = eig(mgamp1);
                mgamp = sort(mgamp);
                mgrid(kcntx,kcnty,:) = [mgamp.'];
            end
        else

            if (~mgridonly)
                [OPTSrlx.v0,erka,flag] = eigs(rka,1,'LM',OPTSrlx);
                OPTSrlx.v0 = real(OPTSrlx.v0); % not sure about this
                if (flag)
                    disp('had trouble converging');
                    amp(kcntx,kcnty) = nan;
                else
                    amp(kcntx,kcnty) = erka;
                end
            end
            
            if (maxlvl > 1)
                [OPTSmg.v0,mgamp,flag] = eigs(mgamp1,1,'LM',OPTSmg);
                OPTSmg.v0 = real(OPTSmg.v0);
                if (flag)
                    disp('had trouble converging');
                    mgrid(kcntx,kcnty) = nan;
                else
                    mgrid(kcntx,kcnty) = mgamp;
                end
                
            end
        end
        
        if (plotfuncs)
            [vmg,emg] = eig(rka);
            [emg,ind] = sort(diag(emg),'descend');
            vmg = vmg(:,ind);

            if (plotfuncs > 1)
                [vmga,emga] = eig(mgamp1);
                [tmp,ind] = sort(diag(emga),'descend');
                vmga = vmga(:,ind);
            end
            
            % PLOT ALL EIGENFUNCTIONS AT THIS WAVENUMBER
            bsz = size(basissave,2);
            syms lx ly
            for vct = 1:howmany
                % SOME TROUBLE WITH SYMBOLICS SO CAN'T USE FUNCTION
                for n = 1:nvar
                    figure
                    evect = vmg(:,vct);
                    plotEigenVect
                    title(['var: ' num2str(n) ' relaxation:' num2str(abs(emg(vct)))]);
                    
                    if (plotfuncs > 1)
                        figure
                        evect = vmga(:,vct);
                        plotEigenVect
                        title(['var: ' num2str(n) ' multigrid:' num2str(abs(emga(vct)))]);
                    end
                end

                % MULTIGRID EFFECT ON RELAXATION EIGENVECTOR
                if (plotfuncs == 3)
                    r1 = create_operator(RB,1,3);
                    k1 = create_operator(KB,1,5);
                    dampedmode = vmg(:,vct) -r1'*mgcycle2d(vw,2,r1*k1*vmg(:,vct));
                    for n = 1:nvar
                       figure
                       evect = dampedmode;
                       plotEigenVect
                       title(['mgrid on fine function: theta/pi:' num2str(kdx/pi)] );
                    end
                end
            end
        end
    end
end

%rlxdamping = max(max(max(abs(amp(:,:,:)-1.0))));
rlxdamping = max(max(max(abs(amp(:,:,:)))));

mgdamping = 0.0;
if (maxlvl > 1)
    mgdamping = max(max(max(abs(mgrid(:,:,:)))));
	[dfactor,IP] = max(mgrid(:,:,size(mgrid,3)));
    if (oned)
        disp('wavenumber in degrees');
        disp([xpos(IP,1)/pi*180]);
    elseif (Fourier2D)
        [dfactor,JP] = max(dfactor);
        IP = IP(JP);
        disp('x and y wavenumber in degrees');
        disp([xpos(IP,JP),ypos(IP,JP)]/pi*180)
    else
        disp('y wavenumber in degrees');
        disp([ypos(1,IP)/pi*180]);
    end
end
    
if (~ploteigs) 
    return
end

figure;
for ic=1:size(amp,3)
    plot(real(amp(:,:,ic)),imag(amp(:,:,ic)),'x');
    hold on;
end
title('relaxation eigenvalues')
axis equal;
hold on
NSID=50;
for cnt=1:NSID
    z(cnt) = +cos(cnt/NSID*2*pi) +sin(cnt/NSID*2*pi)*i;
end
plot(z);

if (maxlvl > 1)
    figure;
    for ic=1:size(mgrid,3)
        plot(real(mgrid(:,:,ic)),imag(mgrid(:,:,ic)),'x');
        hold on
    end
    title('multigrid eigenvalues')
    axis equal;
    hold on
    NSID=50;
    for cnt=1:NSID
        z(cnt) = +cos(cnt/NSID*2*pi) +sin(cnt/NSID*2*pi)*i;
    end
    plot(z);
end


if (ploteigs == 1)
    return;
end

% VERY DETAILED PLOTS
if (~oned && Fourier2D)
    % THIS DOESN'T UNROLL THE 2D EIGENVALUES PROPERLY
    % NOT SURE HOW TO DO IT
    sqrtbsz = sqrt(bsz)
    if (maxlvl > 1)
        figure;
        for ic=1:bsz
            subplot(sqrtbsz,sqrtbsz,ic);
            contourf(xpos,ypos,abs(mgrid(:,:,ic)));
            colorbar;
        end
        title('Multigrid Damping Factor');
    end

    figure;
    for ic=1:bsz
        subplot(sqrtbsz,sqrtbsz,ic);
        contourf(xpos,ypos,abs(amp(:,:,ic)));
        colorbar;
    end
    title('Amplification Factor');
else
    % 1D PLOTTING BEGINS HERE
    % COMPARISON TO ANALYTIC FOR 1D HEAT EQUATION
    if (sys_flag == 0 & rlx_flag==0 && sweep_flag == 0 && omega == -1)
        % figure
        tec_file = ['error_' int2str(P) '_' int2str(scheme) '.dat']
        for ic=1:bsz
            % LOGARITHMIC ERROR PLOTS
            kdxtemp = floor(ic/2)*2*pi*(-1)^(ic+1) +xpos(Nelx/2+1:Nelx,1);
            loglog(abs(kdxtemp),abs(-(amp(Nelx/2+1:Nelx,ic)-1.0)-kdxtemp.^2));
            
            % To save for tecplot
            vector2tecplot(abs(kdxtemp),abs(-(amp(Nelx/2+1:Nelx,ic)-1.0)-kdxtemp.^2),tec_file)
            
            hold on;
            kdxtemp = floor(ic/2)*2*pi*(-1)^(ic+1)+xpos(Nelx/2+1:Nelx,1);
            loglog(abs(kdxtemp),abs(-(amp(Nelx/2+1:Nelx,ic)-1.0)-kdxtemp.^2));
            hold on;
            
            % To save for tecplot
            vector2tecplot(abs(kdxtemp),abs(-(amp(Nelx/2+1:Nelx,ic)-1.0)-kdxtemp.^2),tec_file)
       
        end
        
        xlabel('\theta');
        ylabel('Error in eigenvalue');
        
        figure
        hold on;
        for ic=1:bsz
            % VISUAL COMPARISON TO ANALYTIC
            kdxtemp = floor(ic/2)*2*pi*(-1)^ic +xpos(1:Nelx/2,1);
            plot(kdxtemp,abs(-real(amp(1:Nelx/2,ic)-1.0)));
            plot(kdxtemp,kdxtemp.^2);
            kdxtemp = floor(ic/2)*2*pi*(-1)^(ic+1)+xpos(Nelx/2+1:Nelx,1);
            plot(kdxtemp,abs(-real(amp(Nelx/2+1:Nelx,ic)-1.0)));
            plot(kdxtemp,kdxtemp.^2);
        end
        xlabel('\theta');
        ylabel('eigenvalues');

    elseif (sys_flag == 1 & rlx_flag==0 && sweep_flag == 0 && omega == -1)
        % COMPARISON TO ANALYTIC FOR CONVECTION
        figure
        for ic=1:bsz
            % VISUAL COMPARISON TO ANALYTIC
            % HAD TO SKIP FIRST POINT BECAUSE IT'S NOT GETTING SORTED RIGHT
            kdxtemp = floor(ic/2)*2*pi*(-1)^ic +xpos(2:Nelx/2,1);
            loglog(abs(kdxtemp),abs(amp(2:Nelx/2,ic)-1.0-i*kdxtemp));
            hold on;
            kdxtemp = floor(ic/2)*2*pi*(-1)^(ic+1)+xpos(Nelx/2+1:Nelx,1);
            loglog(abs(kdxtemp),abs(amp(Nelx/2+1:Nelx,ic)-1.0-i*kdxtemp));
        end 
        xlabel('\theta');
        ylabel('error in eigenvalue');
        
        figure
        hold on;
        for ic=1:bsz
            % VISUAL COMPARISON TO ANALYTIC
            % HAD TO SKIP FIRST POINT BECAUSE IT'S NOT GETTING SORTED RIGHT
            kdxtemp = floor(ic/2)*2*pi*(-1)^ic +xpos(2:Nelx/2,1);
            plot(kdxtemp,imag(amp(2:Nelx/2,ic)-1.0));
            plot(kdxtemp,kdxtemp);
            kdxtemp = floor(ic/2)*2*pi*(-1)^(ic+1)+xpos(Nelx/2+1:Nelx,1);
            plot(kdxtemp,imag(amp(Nelx/2+1:Nelx,ic)-1.0));
            rlx = max(abs(amp))
            plot(kdxtemp,kdxtemp);
        end
        xlabel('\theta');
        ylabel('eigenvalues');
    end

    % RELAXATION DAMPING FACTORS
    figure
    hold on;
    for ic=1:bsz
        kdxtemp = floor(ic/2)*2*pi +xpos(1:Nelx,1);
        plot(abs(kdxtemp),abs(amp(1:Nelx,ic)));
    end
    xlabel('|\theta|');
    ylabel('relaxation damping factor');       
   
    % MULTIGRID DAMPING FACTORS
    if (maxlvl > 1)
        figure
        hold on;
        for ic=1:bsz
            kdxtemp = floor(ic/2)*2*pi +xpos(1:Nelx,1);
            plot(abs(kdxtemp),abs(mgrid(1:Nelx,ic)));
        end  
        xlabel('|\theta|');
        ylabel('multigrid damping factor');  
    end
end






