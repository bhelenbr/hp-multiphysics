figure;
set(gcf,'DefaultLineLineWidth',2.0)
hold on;
for n=1:P+1
    ezplot(basis(n),[-1,1])
end
title('Lagrange','FontSize',14);
%title('Spectral Element','FontSize',14);

xlabel('x','Fontsize',14);
h = gca;
set(h,'FontSize', 14,'DataAspectRatioMode','manual','DataAspectRatio',[2 1 1])
axis([-1 1 -.6 1.2]);
return;