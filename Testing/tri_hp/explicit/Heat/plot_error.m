clear;close all
Ngrid = 4:6

strings = {'r','b','g'}
T = 1/(4*pi^2);
for l2p=0:2
    for G=Ngrid
        File = sprintf('Results/L2P_%d_G%d/errs.dat',l2p,G-l2p);
        E = load(File);
        dt = T/(size(E,1)-1);
        semilogy(0:dt:T,E(:,2),strings{l2p+1})
        hold on
    end
end

figure
for l2p=0:2
    File = sprintf('Results/cnvg%d.dat',l2p);
    E = load(File);
    loglog(0.5.^Ngrid,E(:,1),strings{l2p+1})
    hold on
    loglog(0.5.^Ngrid,E(:,2),strings{l2p+1})
    log2(E(2,1)/E(3,1))
    log2(E(2,2)/E(3,2))
end


