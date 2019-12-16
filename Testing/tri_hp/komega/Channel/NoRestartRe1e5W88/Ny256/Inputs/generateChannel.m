% This code generates the geometry file and a restart file for a channel 
% flow driven by a pressure-gradient. This is to be used with
% komega model in tri_hp.
% Run the MATLAB code to solve the channel flow first with the same physical
% constants before running this code

% Physical constants
H = 1; % width of channel section
rho = 1; % density 
ubar_ref = 1; % mean velocity of fully developed channel flow
mu = 1e-5; % dynamic viscosity

% Grid
Nx = 4; % number of elements used in the streamwise direction (y axis) 
Ny = 256; % number of elements used in the cross-stream direction (x axis)
yplus1 = 1; % the  y+ distance of the first node adjacent to the wall
nestedGrid = 1; % 1:finest grid, 2: coarser by a factor of 2, 4: coarser by a factor of 4, ...
wall_res = 0.01;
center_res = 0.01; 

%% Solve for friction factor to find ustar and ystar
Dh = 2*H; % hydraulic diameter
nu = mu/rho; % kinematic viscosity
Re = ubar_ref*Dh/nu; % Reynolds number
syms fr_s;
eq2 = 1/sqrt(fr_s) +2*log10(2.51/(Re*sqrt(fr_s))); % Colebrook's Eq.
fr = double(solve(eq2));
dpdx = -1/2*rho*ubar_ref^2*fr/Dh; % pressure gradient
tau = -dpdx*Dh/4;
ustar = sqrt(tau/rho);
ystar = nu/ustar;

%% Grid Generation with stretching near the wall
xi = (0:1/Ny:1)';

% solving for delta_s which produces the desirable yplus1
stretchingfun = @(delta_s) yplus1*ystar - H/2*(1 +tanh(delta_s*(xi(2)-1/2))/tanh(delta_s/2));
delta = fzero(stretchingfun,5)
y = H/2*(1 +tanh(delta*(xi-1/2))/tanh(delta/2)); % Stretching function (Vinokur 1983)

if nestedGrid~=1
    Ny = Ny/nestedGrid;
    y = y(1:nestedGrid:end);
end

%% Creting geometry and initialization files
% Create .d file
FP = fopen('channel8e7_512.d','w');
fprintf(FP,'%d\n',(Nx+1)*(Ny+1));
count = 0;

% Create initial condition file
FI = fopen('rstrt1_d0_b0.txt','w'); % If using DIRK4 you need to create two files. You can run it twice or ducplicate this one and rename to rstrt1_d1_b0.txt
fprintf(FI,'p0 = 1\n');
fprintf(FI,'npnt = 231, nseg = 614, ntri = 384\n'); % These are wrong but it doesn't matter (?)
fprintf(FI,'END OF HEADER\n');

u_wall = 0; v_wall = 0; ktld_wall = 0; p_wall = 0; % initial conditions for u, v, ktld, p at wall

%% Boundary Points have to be first (counterclock-wise)

% left wall (going from (0,0) to (0,-1))
fprintf(FP,'%d: %21.15e %21.15e %21.15e 1\n',count,0,0,wall_res); % top left corner identified as special to keep pressure constant
count = count+1;
fprintf(FI,'%21.15e %21.15e %21.15e %21.15e %21.15e\n',u_wall, v_wall, ktld_wall, omgtld_wall, p_wall);
for i=1:Nx-1
    fprintf(FP,'%d: %21.15e %21.15e %21.15e 0\n',count,0,-i/Nx,wall_res);
    count = count+1;
    fprintf(FI,'%21.15e %21.15e %21.15e %21.15e %21.15e\n',u_wall, v_wall, ktld_wall, omgtld_wall, p_wall);
end
fprintf(FP,'%d: %21.15e %21.15e %21.15e 2\n',count,0,-1,wall_res);
count = count+1;
fprintf(FI,'%21.15e %21.15e %21.15e %21.15e %21.15e\n',u_wall, v_wall, ktld_wall, omgtld_wall, p_wall);

% bottom periodic boundary: (0,-1) to (1,-1))
for i=2:Ny
    fprintf(FP,'%d: %21.15e %21.15e %21.15e 0\n',count,y(i),-1,center_res);
    count = count+1;
    fprintf(FI,'%21.15e %21.15e %21.15e %21.15e %21.15e\n',0, uc(i), ktldc(i), omgtldc(i), p_wall);
 end

% right wall: (1,-1) to (1,0)
for i = 0:Nx
    fprintf(FP,'%d: %21.15e %21.15e %21.15e 0\n',count,H,-1.0+i/Nx,wall_res);
    count = count +1;
    fprintf(FI,'%21.15e %21.15e %21.15e %21.15e %21.15e\n',u_wall, v_wall, ktld_wall, omgtld_wall, p_wall);
 end

% top periodic boundary: (1,0) to the grid point just before (0,0)
for i=Ny:-1:2
    fprintf(FP,'%d: %21.15e %21.15e %21.15e 0\n',count,y(i),0,center_res);
    count = count+1;
    fprintf(FI,'%21.15e %21.15e %21.15e %21.15e %21.15e\n',0, uc(i), ktldc(i), omgtldc(i), p_wall);
end

%% Now Interior Points
for i=2:Ny
    for j = 1:Nx-1
        fprintf(FP,'%d: %21.15e %21.15e %21.15e 0\n',count,y(i),-j/Nx,center_res);
        count = count+1;
        fprintf(FI,'%21.15e %21.15e %21.15e %21.15e %21.15e\n',0, uc(i), ktldc(i), omgtldc(i), p_wall);
    end
end
fprintf(FI,'b0_s1 plain\nb0_s2 plain\nb0_s1 plain\nb0_s2 plain\nb0_v1 plain\nb0_v2 plain\n');
fclose(FI);

%% Defining the line segments making up the outer boundary
fprintf(FP,'%d\n',2*(Nx+Ny));
count = 0;
%left wall
for i=0:(Nx-1)
    fprintf(FP,'%d: %d %d %d\n',count,i,i+1,1);
    count = count+1;
end

%bottom periodic
for i=Nx:(Ny+Nx-1)
    fprintf(FP,'%d: %d %d %d\n',count,i,i+1,2);
    count = count+1;
end

%right wall
for i=(Nx+Ny):(Ny+2*Nx-1)
    fprintf(FP,'%d: %d %d %d\n',count,i,i+1,1);
    count = count+1;
end

%top periodic
for i=(Ny+2*Nx):(2*Ny+2*Nx-2)
    fprintf(FP,'%d: %d %d %d\n',count,i,i+1,2);
    count = count+1;
end

%last point on top periodic boundary
fprintf(FP,'%d: %d %d %d\n',count,2*Ny+2*Nx-1,0,2);
fclose(FP);





