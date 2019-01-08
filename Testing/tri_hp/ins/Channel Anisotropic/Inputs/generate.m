% Physical constants
H = 0.5; % width of channel section
rho = 1; % density 
ubar = 1; % mean velocity of fully developed channel flow
mu = 1e-4; % dynamic viscosity

% Grid
Ny = 4; % number of elements used in the streamwise direction (y axis) 
Nx = 32; % number of elements used in the cross-stream direction (x axis)
yplus1 = 1; % the  y+ distance of the first node adjacent to the wall
wall_res = 0.01;
center_res = 0.01; 

%% Derived Constants 
Dh = 2*H; % hydraulic diameter
nu = mu/rho; % kinematic viscosity
Re = ubar*Dh/nu; % Reynolds number

% Solve for friction factor 
syms fr_s;
eq2 = 1/sqrt(fr_s) +2*log10(2.51/(Re*sqrt(fr_s))); % Colebrook's Eq.
fr = double(solve(eq2));
dpdx = -1/2*rho*ubar^2*fr/Dh; % pressure gradient
tau = -dpdx*Dh/4;
ustar = sqrt(tau/rho);
ystar = nu/ustar;

%% Grid Generation with stretching near the wall
xi = (0:1/Nx:1)';

% solving for delta_s which produces the desirable yplus1
stretchingfun = @(delta_s) yplus1*ystar - H/2*(1 +tanh(delta_s*(xi(2)-1/2))/tanh(delta_s/2));
delta = fzero(stretchingfun,5)

y = H/2*(1 +tanh(delta*(xi-1/2))/tanh(delta/2)); % Stretching function (Vinokur 1983)

%% Creting geometry and initialization files
% Create .d file
FP = fopen('channel.d','w');
fprintf(FP,'%d\n',(Nx+1)*(Ny+1));
count = 0;

% Create initial condition file
FI = fopen('rstrt1_d0_b0.txt','w');
fprintf(FI,'p0 = 1\n');
fprintf(FI,'npnt = 231, nseg = 614, ntri = 384\n');
fprintf(FI,'END OF HEADER\n');

% Initial Residuals file
FIR = fopen('InitialResiduals.dat','w');

%% Boundary Points have to be first
u_wall = 0; v_wall = 0; ktld_wall = 0; p_wall = 0;
% left wall
fprintf(FP,'%d: 0.0 %f %f 1\n',count,0,wall_res);
count = count+1;
fprintf(FI,'%14.8e %14.8e %14.8e %14.8e %14.8e\n',u_wall, v_wall, ktld_wall, omgtld_wall, p_wall);
fprintf(FIR,'b0 v: %d %10.3e %10.3e %10.3e %10.3e %10.3e\n',count-1,0, R(1), R(2), R(3),0);
for i=1:Ny-1
    fprintf(FP,'%d: 0.0 %f %f 0\n',count,-i/Ny,wall_res);
    count = count+1;
    fprintf(FI,'%14.8e %14.8e %14.8e %14.8e %14.8e\n',u_wall, v_wall, ktld_wall, omgtld_wall, p_wall);
    fprintf(FIR,'b0 v: %d %10.3e %10.3e %10.3e %10.3e %10.3e\n',count-1,0, R(1), R(2), R(3),0);
end
fprintf(FP,'%d: 0.0 %f %f 2\n',count,-1,wall_res);
count = count+1;
fprintf(FI,'%14.8e %14.8e %14.8e %14.8e %14.8e\n',u_wall, v_wall, ktld_wall, omgtld_wall, p_wall);
fprintf(FIR,'b0 v: %d %10.3e %10.3e %10.3e %10.3e %10.3e\n',count-1,0, R(1), R(2), R(3),0);

% bottom periodic boundary
for i=2:Nx
    fprintf(FP,'%d: %f -1.0 %f 0\n',count,y(i),center_res);
    count = count+1;
    fprintf(FI,'%14.8e %14.8e %14.8e %14.8e %14.8e\n',0, uc(i), ktldc(i), omgtldc(i), p_wall);
    fprintf(FIR,'b0 v: %d %10.3e %10.3e %10.3e %10.3e %10.3e\n',count-1,0, R(3*i-2), R(3*i-1), R(3*i),0);
end

% right wall
for i = 0:Ny
    fprintf(FP,'%d: 0.5 %f %f 0\n',count,-1.0+i/Ny,wall_res);
    count = count +1;
    fprintf(FI,'%14.8e %14.8e %14.8e %14.8e %14.8e\n',u_wall, v_wall, ktld_wall, omgtld_wall, p_wall);
    fprintf(FIR,'b0 v: %d %10.3e %10.3e %10.3e %10.3e %10.3e\n',count-1,0, R(3*i+1), R(3*i+2), R(3*i+3),0);
end

% top periodic boundary
for i=Nx:-1:2
    fprintf(FP,'%d: %f 0.0 %f 0\n',count,y(i),center_res);
    count = count+1;
    fprintf(FI,'%14.8e %14.8e %14.8e %14.8e %14.8e\n',0, uc(i), ktldc(i), omgtldc(i), p_wall);
    fprintf(FIR,'b0 v: %d %10.3e %10.3e %10.3e %10.3e %10.3e\n',count-1,0, R(3*i-2), R(3*i-1), R(3*i),0);
end

%% Now Interior Points
for i=2:Nx
    for j = 1:Ny-1
        fprintf(FP,'%d: %f %f %f 0\n',count,y(i),-j/Ny,center_res);
        count = count+1;
        fprintf(FI,'%14.8e %14.8e %14.8e %14.8e %14.8e\n',0, uc(i), ktldc(i), omgtldc(i), p_wall);
        fprintf(FIR,'b0 v: %d %10.3e %10.3e %10.3e %10.3e %10.3e\n',count-1,0, R(3*i-2), R(3*i-1), R(3*i),0);
    end
end
fprintf(FI,'b0_s1 plain\nb0_s2 plain\nb0_s1 plain\nb0_s2 plain\nb0_v1 plain\nb0_v2 plain\n');
fclose(FI);
fclose(FIR);

%% Defining the line segmens making up the outer boundary
fprintf(FP,'%d\n',2*(Ny+Nx));
count = 0;
%left wall
for i=0:(Ny-1)
    fprintf(FP,'%d: %d %d %d\n',count,i,i+1,1);
    count = count+1;
end

%bottom periodic
for i=Ny:(Nx+Ny-1)
    fprintf(FP,'%d: %d %d %d\n',count,i,i+1,2);
    count = count+1;
end

%right wall
for i=(Ny+Nx):(Nx+2*Ny-1)
    fprintf(FP,'%d: %d %d %d\n',count,i,i+1,1);
    count = count+1;
end

%top periodic
for i=(Nx+2*Ny):(2*Nx+2*Ny-2)
    fprintf(FP,'%d: %d %d %d\n',count,i,i+1,2);
    count = count+1;
end

%last point on top periodic boundary
fprintf(FP,'%d: %d %d %d\n',count,2*Nx+2*Ny-1,0,2);
fclose(FP);





