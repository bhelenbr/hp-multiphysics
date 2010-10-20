%function [mgdamping, rlxdamping] = dgdriver(input_vector)
clear global
% To run type: [rlx,mg] = dgdriver([])
% See kdxyloop for changing output options
% FOR DEBUGGING PURPOSES CAN COMMENT ABOVE LINE
% Then can type "clear"
% THEN DEFINE "input_vector = []"
% THEN TYPE "dgdriver"  This leaves all variables in workspace 
% so that they can be queried easily

input_counter = 1;
input_true = size(input_vector,1);


% Polynomial stuff
if (~input_true)
	P = input([num2str(input_counter) ' P:']);
	input_vector(input_counter) = P;
	input_counter = input_counter +1;
else
	P = input_vector(input_counter);
	disp([num2str(input_counter) ' P:' num2str(P)]);
	input_counter = input_counter +1;
end

disp('basis 0 dubiner')
disp('basis 1 legendre')
disp('basis 2 lumped')
disp('basis 3 monomial')
disp('basis 4 nodal')
disp('basis 5 integrated legendre')
if (~input_true)
	basis_flag = input([num2str(input_counter) ' basis:']);
	input_vector(input_counter) = basis_flag;
	input_counter = input_counter +1;
else
	basis_flag = input_vector(input_counter);
	disp([num2str(input_counter) ' basis:' num2str(basis_flag)]);
	input_counter = input_counter +1;
end


% Dimensionality 
global oned;
if (~input_true)
	oned = input([num2str(input_counter) ' 1D? 0/1:']);
	input_vector(input_counter) = oned;
	input_counter = input_counter +1;
else
	oned = input_vector(input_counter);
	disp([num2str(input_counter) ' 1D? 0/1:' num2str(oned)]);
	input_counter = input_counter +1;
end

global Fourier2D;
if (~oned)
	if (~input_true)
		Fourier2D = input([num2str(input_counter) ' 2D Fourier? 0/1:']);
		input_vector(input_counter) = Fourier2D;
		input_counter = input_counter +1;
	else
		Fourier2D = input_vector(input_counter);
		disp([num2str(input_counter) ' 2D Fourier? 0/1:' num2str(Fourier2D)]);
		input_counter = input_counter +1;
	end
end


% SYSTEM INPUTS 
disp('sys_flag = 0 diffusion')
disp('sys_flag = 1 convection')
disp('sys_flag = 2 advection-diffusion')
disp('sys_flag = 3 linearized euler')
if (~input_true)
	sys_flag = input([num2str(input_counter) ' system:']);
	input_vector(input_counter) = sys_flag;
	input_counter = input_counter +1;
else
	sys_flag = input_vector(input_counter);
	disp([num2str(input_counter) ' system:' num2str(sys_flag)]);
	input_counter = input_counter +1;
end

global angle;
angle = 0.0;
if (sys_flag ~= 0 && ~oned) 
	if (~input_true)
		angle = input([num2str(input_counter) ' convection angle:']);
		input_vector(input_counter) = angle;
		input_counter = input_counter +1;
	else
		angle = input_vector(input_counter);
		disp([num2str(input_counter) ' angle:' num2str(angle)]);
		input_counter = input_counter +1;
    end	
    angle = angle/180.0*pi;
end


if (sys_flag == 2) 
	if (~input_true)
		Pe = input([num2str(input_counter) ' Peclet number:']);
		input_vector(input_counter) = Pe;
		input_counter = input_counter +1;
	else
		Pe = input_vector(input_counter);
		disp([num2str(input_counter) ' Pe:' num2str(Pe)]);
		input_counter = input_counter +1;
	end	
end

if (sys_flag == 3) 
	if (~input_true)
		mach = input([num2str(input_counter) ' Mach number:']);
		input_vector(input_counter) = mach;
		input_counter = input_counter +1;
	else
		mach = input_vector(input_counter);
		disp([num2str(input_counter) ' Mach:' num2str(mach)]);
		input_counter = input_counter +1;
	end	
end

if (sys_flag == 0 || sys_flag == 2) 
	% SCHEME LIST
	disp('1 Central-Difference');
	disp('2 Brezzi-et al. 22');
	disp('3 LDG CD');
	disp('4 LDG UPWIND');
	disp('5 Interior Penalty [50]');
	disp('6 Bassi et al. [13]');
	disp('7 Upwind');
	disp('8 Swapping');
	if (~input_true)
		scheme = input([num2str(input_counter) ' scheme?\n']);
		input_vector(input_counter) = scheme;
		input_counter = input_counter +1;
		
		rswp = 0;
		if (scheme == 8)
			rswp = 1;
			scheme = input([num2str(input_counter) ' scheme for swapping? (1-7)']);
			input_vector(input_counter) = scheme;
			input_counter = input_counter +1;	
		end
		
		if (scheme ~= 1 && scheme ~= 7)
			eta = input([num2str(input_counter) ' input \eta/\eta_0:']);
			input_vector(input_counter) = eta;
			input_counter = input_counter +1;
		end
	else
		scheme = input_vector(input_counter);
		disp([num2str(input_counter) ' scheme:' num2str(scheme)]);
		input_counter = input_counter +1;
		
		rswp = 0;
		if (scheme == 8)
			rswp = 1;
			scheme = input_vector(input_counter);
			disp([num2str(input_counter) ' scheme:' num2str(scheme)]);
			input_counter = input_counter +1;
		end
		
		if (scheme ~= 1 && scheme ~= 7)
			eta = input_vector(input_counter);
			disp([num2str(input_counter) ' input \eta/\eta_0:' num2str(eta)]);
			input_counter = input_counter +1;
		end
		
	end		
end

prcndtn = 0;
if (sys_flag == 3)
	if (~input_true)
		prcndtn = input([num2str(input_counter) ' use Turkel preconditioner to normalize wave speeds?:']);
		input_vector(input_counter) = prcndtn;
		input_counter = input_counter +1;
	else
		prcndtn = input_vector(input_counter);
		disp([num2str(input_counter) ' prcndtn:' num2str(prcndtn)]);
		input_counter = input_counter +1;
	end	
end

ksupg = 0;
if (sys_flag ~= 0)
	if (~input_true)
		ksupg = input([num2str(input_counter) ' add supg upwinding terms to stiffness matrix?:']);
		input_vector(input_counter) = ksupg;
		input_counter = input_counter +1;
	else
		ksupg = input_vector(input_counter);
		disp([num2str(input_counter) ' ksupg:' num2str(ksupg)]);
		input_counter = input_counter +1;
	end	
end







% RELAXATION SCHEME STUFF
global rlx_flag
disp('rlx_flag = 0 mass matrix')
disp('rlx_flag = 1 jacobi')
disp('rlx_flag = 2 element jacobi')
disp('rlx_flag = 3 jacobi of m^-1 k')
disp('rlx_flag = 4 Line Solve');
disp('rlx_flag = 5 Alternating Direction Line Solve (2D only)');
disp('rlx_flag = 6 Alternating Direction Implicit (2D only)');
disp('rlx_flag = 7 incomplete LU (only 2D)')
if (~input_true)
	rlx_flag = input([num2str(input_counter) ' relaxation scheme:']);
	input_vector(input_counter) = rlx_flag;
	input_counter = input_counter +1;
else
	rlx_flag = input_vector(input_counter);
	disp([num2str(input_counter) ' relaxation scheme:' num2str(rlx_flag)]);
	input_counter = input_counter +1;
end

sweep_flag = 0;
if (rlx_flag < 6)
	disp('sweep_flag = 0 no sweeping');
	disp('sweep_flag = 1 Block Gauss-Seidel');
	disp('sweep_flag = 2 Symmetric Block Gauss Seidel');
	if (~input_true)
		sweep_flag = input([num2str(input_counter) ' sweeping:']);
		input_vector(input_counter) = sweep_flag;
		input_counter = input_counter +1;
	else
		sweep_flag = input_vector(input_counter);
		disp([num2str(input_counter) ' sweeping:' num2str(sweep_flag)]);
		input_counter = input_counter +1;
	end
end

if (rlx_flag == 6 && sweep_flag == 2)
    disp('SGS ADL not implemented');
end

% ANALYSIS TYPES: 0 = 1-SWEEP, 1 = 2-SWEEP RESIDUAL EVALUATION, 2 = 2-SWEEP NO RESIDUAL
global analysis_type
analysis_type = 0;
if (sweep_flag == 2) 
    analysis_type = 1;
end
if (rlx_flag == 5) 
    analysis_type = 1;
end
if (rlx_flag == 6 || rlx_flag == 7)
    analysis_type = 2;
end


implicit_flag = 0;
if (rlx_flag ~= 6 && rlx_flag ~= 0)
	if (~input_true)
		implicit_flag = input([num2str(input_counter) ' add implicit term?:']);
		input_vector(input_counter) = implicit_flag;
		input_counter = input_counter +1;
	else
		implicit_flag = input_vector(input_counter);
		disp([num2str(input_counter) ' add implicit term?:' num2str(implicit_flag)]);
		input_counter = input_counter +1;
	end	
end

if (implicit_flag || rlx_flag==6) 
	if (~input_true)
		 dt = input([num2str(input_counter) ' time step?:']);
		input_vector(input_counter) = dt;
		input_counter = input_counter +1;

		mu = input([num2str(input_counter) ' mu implicit factor (0.0 - 1.0)?:']);
		input_vector(input_counter) = mu;
		input_counter = input_counter +1;		
		
	else
		dt = input_vector(input_counter);
		disp([num2str(input_counter) ' time step?:' num2str(dt)]);
		input_counter = input_counter +1;
		
		mu = input_vector(input_counter);
		disp([num2str(input_counter) ' mu implicit factor (0.0 - 1.0)?:' num2str(mu)]);
		input_counter = input_counter +1;
	end	
end

msupg = 0;
if (ksupg && (implicit_flag || rlx_flag == 0))
	if (~input_true)
		msupg = input([num2str(input_counter) ' add supg upwinding terms to mass matrix?:']);
		input_vector(input_counter) = msupg;
		input_counter = input_counter +1;
	else
		msupg = input_vector(input_counter);
		disp([num2str(input_counter) ' msupg:' msupg]);
		input_counter = input_counter +1;
	end	
end
global omega
omega = -1;
if (~implicit_flag) 
	if (~input_true)
		omega = input([num2str(input_counter) ' overrelaxation factor?\n (less than zero will turn off scaling by maximum eigenvalue)']);
		input_vector(input_counter) = omega;
		input_counter = input_counter +1;
	else
		omega = input_vector(input_counter);
		disp([num2str(input_counter) ' omega:' num2str(omega)]);
		input_counter = input_counter +1;
	end	
end

global rk_flag;
rk_flag = 0;
if (~input_true)
	rk_flag = input([num2str(input_counter) ' use Runge Kutta Scheme?']);
	input_vector(input_counter) = rk_flag;
	input_counter = input_counter +1;
else
	rk_flag = input_vector(input_counter);
	disp([num2str(input_counter) ' use Runge Kutta Scheme?' num2str(rk_flag)]);
	input_counter = input_counter +1;
end	


% GRID STUFF
global Nelx
if (~input_true)
	Nel = input([num2str(input_counter) ' Number of elements']);
	input_vector(input_counter) = Nel;
	input_counter = input_counter +1;
else
	Nel = input_vector(input_counter);
	disp([num2str(input_counter) ' Number of elements' num2str(Nel)]);
	input_counter = input_counter +1;
end	

if (oned) 
    Nelx = Nel;
    Nely = 1;
else
    Nelx = Nel;
    Nely = Nel;
end

dyodx = 1.0;
if (~oned)
	if (~input_true)
		dyodx = input([num2str(input_counter) ' Dy/Dx']);
		input_vector(input_counter) = dyodx;
		input_counter = input_counter +1;
	else
		dyodx = input_vector(input_counter);
		disp([num2str(input_counter) ' Dy/Dx' num2str(dyodx)]);
		input_counter = input_counter +1;
	end	
end
	
% MULTIGRID STUFF
global maxlvl;
if (~input_true)
	maxlvl = input([num2str(input_counter) ' multigrid levels (1 - no multigrid)']);
	input_vector(input_counter) = maxlvl;
	input_counter = input_counter +1;
else
	maxlvl = input_vector(input_counter);
	disp([num2str(input_counter) ' multigrid levels (1 - no multigrid)' num2str(maxlvl)]);
	input_counter = input_counter +1;
end	

cgswitch = 0;
vw = 0;
krestrict = 0;
if (maxlvl > 1)
	if (~input_true)
		krestrict = input([num2str(input_counter) ' restrict to create coarse grid operators 0/1?']);
		input_vector(input_counter) = krestrict;
		input_counter = input_counter +1;
	else
		krestrict = input_vector(input_counter);
		disp([num2str(input_counter) ' restrict to create coarse grid operators 0/1?' num2str(krestrict)]);
		input_counter = input_counter +1;
	end		
	
	if (~input_true)
		vw = input([num2str(input_counter) ' v/w cycle 1/2?']);
		input_vector(input_counter) = vw;
		input_counter = input_counter +1;
	else
		vw = input_vector(input_counter);
		disp([num2str(input_counter) ' v/w cycle 1/2?' num2str(vw)]);
		input_counter = input_counter +1;
    end		
	
    if (P/2^(maxlvl-1) < 1 && P > 0)
        if (~input_true)
            cgswitch = input([num2str(input_counter) ' switch to continuous space?']);
            input_vector(input_counter) = cgswitch;
            input_counter = input_counter +1;
        else
            cgswitch = input_vector(input_counter);
            disp([num2str(input_counter) ' switch to continuous space?' num2str(cgswitch)]);
            input_counter = input_counter +1;
        end	
    end
end

supg_restrict = 0;
if (cgswitch && sys_flag ~= 0)
    if (~input_true)
        supg_restrict = input([num2str(input_counter) ' use supg restriction?']);
        input_vector(input_counter) = supg_restrict;
        input_counter = input_counter +1;
    else
        supg_restrict = input_vector(input_counter);
        disp([num2str(input_counter) ' use supg restriction?' num2str(supg_restrict)]);
        input_counter = input_counter +1;
    end	
end


if ~input_true
    disp('This is your input vector');
    ivect = '[';
    for indx=1:size(input_vector,2)-1
        ivect = [ivect sprintf('%g,',input_vector(indx))];
    end
    ivect = [ivect sprintf('%g]',input_vector(size(input_vector,2)))];
    disp(ivect)
end


global KB MB MB1 RB;
global D Ly Uy Lx Ux;  % FOR ILU

tic
% PRECALCULATE MATRICES FOR EACH LEVEL
Psave = P;
PC = floor(P/2);
%PC = max(P-1,0) %% TEMPORARY
%PC = floor((P+1)/2)-1;  %% TEMPORARY
for lvl = 1:maxlvl
    % CALCULATE STIFFNESS AND RELAXATION MATRICES
    dgstiffness2d;
    P = PC;
    PC = floor(P/2);
    %PC = max(P-1,0) %% TEMPORARY
    %PC = floor((P+1)/2)-1;  %% TEMPORARY
end
toc

P=Psave;
tic
kdxyloop;
toc
