% Michelle Pede 
% Research under Valentine summer 08
% movement of creepy crawler
clc
clear


ke=1.3;
h=0.7;
w=pi;
k=1.3;
dx=1;
len=7;
x0=0:dx:len;
n=length(x0);
dt=.1;
time=0;
t=0:dt:time;
timesteps=length(t);
s=0.1;

x1=zeros(1,n);

deri=zeros(1,n);
thet=zeros(1,n);
xtop=zeros(1,n);
ytop=zeros(1,n);
xbot=zeros(1,n);
ybot=zeros(1,n);

r=1;
c0=-sqrt(r^2-s^2);
y0=0;
m=10; %Don't make this too big or circle will overlap flagella
[xx,yy]=circle(r,c0,y0,m);
ex0=len;
denom=.7;
for i = 1:timesteps
    x1=h*(1-exp(-ke^2*x0.^2)).*cos(k*x0-w*t(i));
    deri=2*h*ke^2*x0.*exp(-ke^2*x0.^2).*cos(k*x0-w*t(i))-h*(1-exp(-ke^2*x0.^2)).*sin(k*x0-w*t(i))*k;
    thet=atan(deri(1,:));
%     ytop2=h*(1-exp(-ke^2*x0.^2)).*cos(k*x0-w*t(i))+s*(1-exp(-(x0-ex0).^2/e^2)).*(1+(2*h*ke^2*x0.*exp(-ke^2*x0.^2).*cos(k*x0-w*t(i))-h*(1-exp(-ke^2*x0.^2)).*sin(k*x0-w*t(i))*k).^2).^(1/2);
    ytop=x1+s*(1-exp(-((x0-ex0)/denom).^2))./cos(thet);
    ybot=x1-s*(1-exp(-((x0-ex0)/denom).^2))./cos(thet);
    
    drawnow
    plot(x0,x1,'k',x0,ytop,'b',x0,ybot,'b',xx,yy,'b');
    axis image
end


file=fopen('flagellum.d','wt');
fprintf(file,'%g \n',length(x0)*2+length(xx)-2+3);

% Outer Box
xmin=-10;
xmax=10;
ymin=-10;
ymax=10;
bsize=1;
fsize=.1;


fprintf(file,'0: %15.15f %15.15f %15.15f 0 \n',xmin,ymin,bsize);
fprintf(file,'1: %15.15f %15.15f %15.15f 0 \n',xmax,ymin,bsize);
fprintf(file,'2: %15.15f %15.15f %15.15f 0 \n',xmax,ymax,bsize);
fprintf(file,'3: %15.15f %15.15f %15.15f 0 \n',xmin,ymax,bsize);
ind = 4;

%top line 
for i=1:length(x0)-1
    fprintf(file,'%g: %15.15f %15.15f %15.15f 0 \n',ind,x0(i),ytop(i),fsize);
    ind=ind+1;
end 
fprintf(file,'%g: %15.15f %15.15f %15.15f %g \n',ind,x0(end),ytop(end),fsize,ind);
endvrtx=ind;
ind=ind+1;
%bottom line
for i=length(x0)-1:-1:1
    fprintf(file,'%g: %15.15f %15.15f %15.15f 0 \n',ind,x0(i),ybot(i),fsize);
    ind=ind+1;
end 
%cell circle
for i=length(xx)-1:-1:2
    fprintf(file,'%g: %15.15f %15.15f %15.15f 0 \n',ind,xx(i),yy(i),fsize);
    ind=ind+1;
end

fprintf(file,'%g \n',length(x0)*2+length(xx)-2+3);

%Connect Outer Box
for i=1:3
fprintf(file,'%g: %g %g %g \n',i-1, i-1, i, i);
end
fprintf(file,'%g: %g %g %g \n',3,3,0,4);

ind=4;
%Connect Top Line
for i=1:length(xtop)-1
    fprintf(file,'%g: %g %g %g \n',ind,ind,ind+1,6);
    ind=ind+1;
end

%Connect bottom
for i=1:length(xbot)-1
    fprintf(file,'%g: %g %g %g \n',ind,ind,ind+1,7);
    ind=ind+1;
end

%Connect end circle
for i=1:length(xx)-2
    fprintf(file,'%g: %g %g %g \n',ind,ind,ind+1,5);
    ind=ind+1;
end
fprintf(file,'%g: %g %g %g \n',ind,ind,4,5);
fclose(file);


% b0_v0_type: symbolic
% b0_v0_locx0:  unsteady with t as time variable
% b0_v0_locx1:  unsteady with t as time variable
% b0_v0_r_type: moving
% 
% b0_v1_type: symbolic
% b0_v1_locx0:  unsteady with t as time variable
% b0_v1_locx1:  unsteady with t as time variable
% b0_v1_r_type: moving
% 
% Also for the flagellum sides, you will need:
% b0_s?_r_type: deforming

file=fopen('createmesh.inpt','wt');

fprintf(file,'b0_s1_type: plain\n');
fprintf(file,'b0_s2_type: plain\n');
fprintf(file,'b0_s3_type: plain\n');
fprintf(file,'b0_s4_type: plain\n');
fprintf(file,'s: %15.15f\n',s);
fprintf(file,'h: %15.15f\n',h);
fprintf(file,'w: %15.15f\n',w);
fprintf(file,'ke: %15.15f\n',ke);
fprintf(file,'k: %15.15f\n',k);
fprintf(file,'denom: %15.15f\n',denom);

fprintf(file,'ex0: %15.15f\n',len);

fprintf(file,'b0_v%g_r_type: moving\n',endvrtx);
fprintf(file,'b0_v%g_type: symbolic\n',endvrtx);
fprintf(file,'b0_v%g_locx0: ex0 \n',endvrtx);
fprintf(file,'b0_v%g_locx1: h*(1-exp(-ke^2*ex0^2))*cos(k*ex0-w*t) \n',endvrtx);
 
%left circle
fprintf(file,'b0_s5_type: circle\n');
fprintf(file,'b0_s5_radius: %15.15f\n',r);
fprintf(file,'c0: -sqrt(b0_s5_radius^2-s^2)\n');
fprintf(file,'b0_s5_center: c0 0\n');


%top of flagellum tail
fprintf(file,'b0_s6_r_type: deforming\n');
fprintf(file,'b0_s6_type: symbolic\n');
fprintf(file,'b0_s6_h: h*(1-exp(-ke^2*x0^2))*cos(k*x0-w*t)+s*(1-exp(-(x0-ex0)^2/denom^2))*(1+(2*h*ke^2*x0*exp(-ke^2*x0^2)*cos(k*x0-w*t)-h*(1-exp(-ke^2*x0^2))*sin(k*x0-w*t)*k)^2)^(1/2)-x1 \n');
fprintf(file,'b0_s6_dhdx0: 2*h*ke^2*x0*exp(-ke^2*x0^2)*cos(k*x0-w*t)-h*(1-exp(-ke^2*x0^2))*sin(k*x0-w*t)*k+2*s*(x0-ex0)/denom^2*exp(-(x0-ex0)^2/denom^2)*(1+(2*h*ke^2*x0*exp(-ke^2*x0^2)*cos(k*x0-w*t)-h*(1-exp(-ke^2*x0^2))*sin(k*x0-w*t)*k)^2)^(1/2)+s*(1-exp(-(x0-ex0)^2/denom^2))/(1+(2*h*ke^2*x0*exp(-ke^2*x0^2)*cos(k*x0-w*t)-h*(1-exp(-ke^2*x0^2))*sin(k*x0-w*t)*k)^2)^(1/2)*(2*h*ke^2*x0*exp(-ke^2*x0^2)*cos(k*x0-w*t)-h*(1-exp(-ke^2*x0^2))*sin(k*x0-w*t)*k)*(2*h*ke^2*exp(-ke^2*x0^2)*cos(k*x0-w*t)-4*h*ke^4*x0^2*exp(-ke^2*x0^2)*cos(k*x0-w*t)-4*h*ke^2*x0*exp(-ke^2*x0^2)*sin(k*x0-w*t)*k-h*(1-exp(-ke^2*x0^2))*cos(k*x0-w*t)*k^2)\n');
fprintf(file,'b0_s6_dhdx1: -1\n');

%bottom of flagellum tail
fprintf(file,'b0_s7_r_type: deforming\n');
fprintf(file,'b0_s7_type: symbolic\n');
fprintf(file,'b0_s7_h: h*(1-exp(-ke^2*x0^2))*cos(k*x0-w*t)-s*(1-exp(-(x0-ex0)^2/denom^2))*(1+(2*h*ke^2*x0*exp(-ke^2*x0^2)*cos(k*x0-w*t)-h*(1-exp(-ke^2*x0^2))*sin(k*x0-w*t)*k)^2)^(1/2)-x1 \n');
fprintf(file,'b0_s7_dhdx0: 2*h*ke^2*x0*exp(-ke^2*x0^2)*cos(k*x0-w*t)-h*(1-exp(-ke^2*x0^2))*sin(k*x0-w*t)*k-2*s*(x0-ex0)/denom^2*exp(-(x0-ex0)^2/denom^2)*(1+(2*h*ke^2*x0*exp(-ke^2*x0^2)*cos(k*x0-w*t)-h*(1-exp(-ke^2*x0^2))*sin(k*x0-w*t)*k)^2)^(1/2)-s*(1-exp(-(x0-ex0)^2/denom^2))/(1+(2*h*ke^2*x0*exp(-ke^2*x0^2)*cos(k*x0-w*t)-h*(1-exp(-ke^2*x0^2))*sin(k*x0-w*t)*k)^2)^(1/2)*(2*h*ke^2*x0*exp(-ke^2*x0^2)*cos(k*x0-w*t)-h*(1-exp(-ke^2*x0^2))*sin(k*x0-w*t)*k)*(2*h*ke^2*exp(-ke^2*x0^2)*cos(k*x0-w*t)-4*h*ke^4*x0^2*exp(-ke^2*x0^2)*cos(k*x0-w*t)-4*h*ke^2*x0*exp(-ke^2*x0^2)*sin(k*x0-w*t)*k-h*(1-exp(-ke^2*x0^2))*cos(k*x0-w*t)*k^2)\n');
fprintf(file,'b0_s7_dhdx1: -1\n');


fprintf(file,'nblock: 1\n');
fprintf(file,'b0_filetype: 8\n');
fprintf(file,'b0_growth factor: 1000.0\n');
fprintf(file,'b0_mesh: ./flagellum\n');
fprintf(file,'b0_adapt: 1\n');

fprintf(file,'#logfile: translate\n');
fprintf(file,'# number of geometric multigrid levels: \n');
fprintf(file,'#ngrid: 3\n');
fprintf(file,'# iterations on coarsening in multigrid cycle: \n');
fprintf(file,'#itercrsn: 3\n');
fprintf(file,'# Number of time steps in simulation: \n');
fprintf(file,'ntstep: 10\n');
fprintf(file,'# number of iterative multigrid cycles per time step: \n');
fprintf(file,'ncycle: 4\n');
fprintf(file,'# inverse time step:\n');
fprintf(file,'dtinv: 20.0\n');
fprintf(file,'out_intrvl: 1\n');

fclose(file);