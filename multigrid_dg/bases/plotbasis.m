figure;
cd bases
P = 4
PC = 2
monomial
set(gcf,'DefaultLineLineWidth',2.0)
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1) pos(2) 0.5*pos(3) pos(4)]);
hold on;
for n=1:P+1
    ezplot(basis(n),[-1,1])
end
title('monomial','FontSize',14);
%title('Gauss-Lobatto-Lagrange','FontSize',14);

xlabel('x','Fontsize',14);
h = gca;
set(h,'FontSize', 14,'DataAspectRatio',[2 1 1],'DataAspectRatioMode','manual')
axis([-1 1 -1 1]);
cd ..;
