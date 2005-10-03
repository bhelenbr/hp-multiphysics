load mgrid2
mgrid2 = mgrid;
load mgrid4
mgrid4 = mgrid;
load mgrid8
mgrid8 = mgrid;
load mgrid16;
mgrid16 = mgrid;

figure;
set(gcf,'DefaultLineLineWidth',2.0)
%set(gcf,'DefaultMarkerSize',2.0)
subplot(2,2,1)
hold on;
P=2;
for i=1:size(mgrid2,2)-1
	plot(mgrid2(:,1),abs(mgrid2(:,i+1)),'x');
end
xlabel('Theta','Fontsize',14);
title('P=2','FontSize',14);
h = gca;
set(h,'FontSize', 14,'DataAspectRatioMode','manual','DataAspectRatio',[1 1/2 1])
axis([0 pi 0 1]);

subplot(2,2,2)
hold on;
P=4;
for i=1:size(mgrid4,2)-1
	plot(mgrid4(:,1),abs(mgrid4(:,i+1)),'x');
end
xlabel('Theta','Fontsize',14);
title('P=4','FontSize',14);
h = gca;
set(h,'FontSize', 14,'DataAspectRatioMode','manual','DataAspectRatio',[1 1/2 1])
axis([0 pi 0 1]);

subplot(2,2,3)
hold on;
P=8;
for i=1:size(mgrid8,2)-1
	plot(mgrid8(:,1),abs(mgrid8(:,i+1)),'x');
end
xlabel('Theta','Fontsize',14);
title('P=8','FontSize',14);
h = gca;
set(h,'FontSize', 14,'DataAspectRatioMode','manual','DataAspectRatio',[1 1/2 1])
axis([0 pi 0 1]);

subplot(2,2,4)
hold on;
P=16;
for i=1:size(mgrid16,2)-1
	plot(mgrid16(:,1),abs(mgrid16(:,i+1)),'x');
end
xlabel('Theta','Fontsize',14);
title('P=16','FontSize',14);
h = gca;
set(h,'FontSize', 14,'DataAspectRatioMode','manual','DataAspectRatio',[1 1/2 1])
axis([0 pi 0 1]);