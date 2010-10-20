clear
%close all

% whattodo = 1: only angle study or Mach study (1 variable loop)
% whattodo = 2: ILU dt study (2 variable loop)
% whattodo = 3: grid & p, and angle study (3 variable loop)
% whattodo = 4: (4 variable loop)

whattodo = 1

fname = input('filename: (put in single quotes) ');


if (whattodo == 1)
    p = 0;
    bas = 3;
    oned = 0;
    Fr2d = 0;
    sys = 3;
    thetas = 36;
    mach = 0.1;
    prcndtn = 0;
    addsupg = 0;
    
    % BJ (rlx = 2) (SGS: swp = 2)
    % LINE SOLVE (rlx = 4) (SGS: swp = 2)
    % ADL (rlx = 5) GS swp = 1
    % ILU (rlx = 7)        
    rlxscheme = 5;
    swp = 1; 
    implct = 0; 
    dts = 10.0;
    mu = 1.0;
    omega = -2/3;
    
    rk = 0;
    grd = 8;
    aratio = 1;
    mglvl = 2;
    rstrct = 1;
    vw = 1;
    cgswitch = 0;
    
     thetas = [0:15:360];
%     thetas = [0,36,135];
     nstudy = length(thetas);
%    mach = [0.01,0.05, 0.1, 0.5];
%    nstudy = length(mach);
%     dts = [30:10:50];
%     nstudy = length(dts);
%    grd = [8,16,32];
%    nstudy = length(grd);

    mg = zeros(length(nstudy));
    rlx = zeros(length(nstudy));
    for tc=1:nstudy
        thc = min(length(thetas),tc);
        mc = min(length(mach),tc);
        dtc = min(length(dts),tc);
        gc = min(length(grd),tc);

        if (implct) 
            [mg(tc),rlx(tc)] = dgdriver([p,bas,oned,Fr2d,sys,thetas(thc),mach(mc),prcndtn,addsupg,rlxscheme,swp,implct,dts(dtc),mu,rk,grd(gc),aratio,mglvl,rstrct,vw,cgswitch])
        elseif (rlxscheme ~= 7)
            [mg(tc),rlx(tc)] = dgdriver([p,bas,oned,Fr2d,sys,thetas(thc),mach(mc),prcndtn,addsupg,rlxscheme,swp,implct,omega,rk,grd(gc),aratio,mglvl,rstrct,vw,cgswitch])
        else
            [mg(tc),rlx(tc)] = dgdriver([p,bas,oned,Fr2d,sys,thetas(thc),mach(mc),prcndtn,addsupg,rlxscheme,implct,omega,rk,grd(gc),aratio,mglvl,rstrct,vw,cgswitch])
        end
        
        tc/nstudy
    end


    
    plot(thetas,mg, 'k-')
    hold on
    plot(thetas,rlx, 'k--')
    xlabel('Flow Angle/Degrees');
    ylabel('Damping Factor');
    axis([0 360, 0.25 1.05])
    
%    semilogx(mach,mg,'k-')
%    xlabel('Mach Number');
%    ylabel('DF','Rotation',0);

%     plot(dts,rlx,'k-')
%     hold on
%     plot(dts,mg,'b-')
%     xlabel('Time step,\Deltat','FontSize', 14, 'FontName', 'Times');
%     ylabel('Damping Factor','FontSize', 14,'FontName','Times');


%      semilogx(grd,mg,'kx-')
%      hold on
%      xlabel('N');
%      ylabel('Damping Factor');
%      axis([grd(1) grd(tc) 0.25 1.0]);

    % save ILU_dt_study
    % save line_implicit_theta36
    % save line_w23_theta36
    
    fhandle = gcf;
    saveas(fhandle,['../Figures/' fname '.fig']);
    system(['cp looper.m ../Figures/' fname '.m']);
    
    
    return
end

if (whattodo == 2)
    ps = [1,2,4];
    %ps = [2];
    bas = 3;
    oned = 0;
    Fr2d = 0;
    sys = 3;
    %thetas = [0:5:50];
    thetas = [0,36,144,180];
    mach = 0.1;
    prcndtn = 0.0;
    addsupg = 0;
    rlxscheme = 7;
    swp = 0;
    implct = 0;
    omega = -0.5;
    rk = 0;
    grds = [8,16];
    aratio = 1;
    mglvl = 2;
    rstrct = 1;
    vw = 1;
    cgswitch = 0;

    mg = zeros(length(ps),length(grds),length(thetas));
    rlx = zeros(length(ps),length(grds),length(thetas));

    for pc=1:size(ps,2)
        for gc=1:size(grds,2)
            for tc=1:size(thetas,2)
                % ILU STUDY
                 [mg(pc,gc,tc),rlx(pc,gc,tc)] = dgdriver([ps(pc),bas,oned,Fr2d,sys,thetas(tc),mach,prcndtn,addsupg,rlxscheme,swp,omega,rk,grds(gc),aratio,mglvl,rstrct,vw,cgswitch])
            end
        end
    end

    save ILU_w.5
    return
end

if (whattodo == 3)
    %ps = [1,2,4];
    ps = [2];
    bas = 3;
    oned = 0;
    Fr2d = 0;
    sys = 3;
    %thetas = [0:5:50];
    thetas = [0,36,144,180];
    thetas = [-15,0,15,30,45,75,85,95,105,120,135,150,165,180,195];
    thetas = [144];
    mach = 0.1;
    prcndtn = 0.0;
    addsupg = 0;
    rlxscheme = 7;
    swp = 0;
    implct = 0;
    dt = 3/mach;
    mu = 1.0
    omega = -0.5;
    rk = 0;
    grds = [8,16,32];
    aratio = 1;
    mglvl = 2;
    rstrct = 1;
    vw = 1;
    cgswitch = 0;

    mg = zeros(length(thetas),length(grds),length(ps));
    rlx = zeros(length(thetas),length(grds),length(ps));

    for pc=1:length(ps)
        for gc=1:length(grds)
            for tc=1:length(thetas)
                % ILU STUDY Under-relaxed
                [mg(tc,gc,pc),rlx(tc,gc,pc)] = dgdriver([ps(pc),bas,oned,Fr2d,sys,thetas(tc),mach,prcndtn,addsupg,rlxscheme,implct,omega,rk,grds(gc),aratio,mglvl,rstrct,vw,cgswitch])
                % ILU study ITA
                % [mg(tc,gc,pc),rlx(tc,gc,pc)] = dgdriver([ps(pc),bas,oned,Fr2d,sys,thetas(tc),mach,prcndtn,addsupg,rlxscheme,implct,dt,mu,rk,grds(gc),aratio,mglvl,rstrct,vw,cgswitch]);
            end
        end
    end
    
    % To plot against third variable need to do something like this 
    % *stupid matlab thing 
    % mg = permute(mg,[3,1,2]);
    % rlx = permute(rlx,[3,1,2]);
%     for pc=1:length(ps)
%         for gc=1:length(grds)
%             plot(thetas,mg(:,gc,pc), 'k-')
%             hold on
%             plot(thetas,rlx(:,gc,pc), 'k--')
%             xlabel('Flow Angle/Degrees');
%             ylabel('DF','Rotation',0);  % Stupid matlab bug!!
%             axis([thetas(1) thetas(length(thetas)) 0.25 1.05])
%             text(thetas(2),mg(2,gc,pc),['p=' int2str(ps(pc))])
%         end
%     end
    
    for pc=1:length(ps)
        for tc=1:length(thetas)
            semilogx(grds,mg(tc,:,pc), 'k-')
            hold on
            semilogx(grds,rlx(tc,:,pc), 'k--')
            xlabel('Grid Resolution');
            ylabel('DF','Rotation',0);  % Stupid matlab bug!!
            axis([grds(1) grds(length(grds)) 0.25 1.05])
            text(grds(2),mg(tc,2,pc),['p=' int2str(ps(pc))])
        end
    end

    return
end

if (whattodo == 4)
    ps = [1,2,4];
    ps = [2];
    bas = 3;
    oned = 0;
    Fr2d = 0;
    sys = 3;
    thetas = [0:5:50];
    thetas = [0,36];
    thetas = [0];
    mach = 0.1;
    prcndtn = 0;
    addsupg = 0;
    rlxscheme = 2;
    swp = 0;
    implct = 1;
    mu = 1.0;
    dts =  [1. 4. 8.0 12.0 16.0 24.0 32.0 40.0 48.0 64.0];
    dts = [100];
    omega = -1;
    rk = 0;
    grds = [16,32];
    aratio = 1;
    mglvl = 2;
    rstrct = 1;
    vw = 1;
    cgswitch = 0;
    
    for pc=1:size(ps,2)
        for gc=1:size(grds,2)
            for tc=1:size(thetas,2)
                for dtc=1:size(dts,2);
                    % BJ STUDY
                    [mg(pc,gc,tc,dtc),rlx(pc,gc,tc,dtc)] = dgdriver([ps(pc),bas,oned,Fr2d,sys,thetas(tc),mach,prcndtn,addsupg,rlxscheme,swp,implct,dts(dtc),mu,rk,grds(gc)/ps(pc),aratio,mglvl,rstrct,vw,cgswitch])
                    cycperdec(dtc) = log(0.1)/log(mg(pc,gc,tc,dtc))
                end
                figure
                semilogx(dts,cycperdec);
                atitle = ['p=' num2str(ps(pc)) ' Nel=' num2str(grds(gc)) ' theta=' num2str(thetas(tc))];
                myexport(atitle);
            end
        end
    end
end
