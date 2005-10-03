clear all;
disp('input 1: P');
disp(' ');

global sys_flag;
disp('input 2: sys_flag?');
disp('sys_flag = 0 convection')
disp('sys_flag = 1 diffusion')
disp('sys_flag = 2 linearized euler')
disp(' ');

disp('input 3: basis_flag?');
disp('basis_flag = 0 dubiner')
disp('basis_flag = 1 legendre')
disp('basis_flag = 2 lumped')
disp('basis_flag = 3 monomial')
disp('basis_flag = 4 nodal')
disp(' ');

global rlx_flag;
disp('input 4: relaxation scheme?');
disp('rlx_flag = 0 mass matrix')
disp('rlx_flag = 1 jacobi')
disp('rlx_flag = 2 element jacobi')
disp('rlx_flag = 3 jacobi of m^-1 k')
disp(' ');

disp('input 5: sweep_flag = 0/1');
disp(' ');

global omega;
disp('input 6: overrelaxation factor?');
disp(' ');

disp('input 7: restrict to create coarse grid operators 0/1?');
disp(' ');

global maxlvl;
disp('input 8: multigrid levels (1 - no multigrid)');
disp(' ');

disp('input 9: v/w cycle 1/2?');
disp(' ');

disp('input 10: onedimensional? 0/1');
disp(' ');

vinput = input('give vector for inputs in form [....]');
P = vinput(1)
Pfine = P;
sys_flag = vinput(2)
basis_flag = vinput(3)
rlx_flag = vinput(4)
sweep_flag = vinput(5)
omega = vinput(6)
krestrict = vinput(7)
if ( P > 0 )
    maxlvl = min(vinput(8),log2(P)+2)
else
    maxlvl = 1
end
vw = vinput(9)
oned = vinput(10)


if (sys_flag ~= 1)
	ksupg = input('add supg upwinding terms to stiffness?\n');
	msupg = input('add supg upwinding terms to mass matrix?\n');
end

global KB MB RB vrestrict vmaxeig nvar;

% PRECALCULATE MATRICES FOR EACH LEVEL
PC = floor(P/2);
count = 1;
for lvl = 1:maxlvl
    
    % CALCULATE STIFFNESS AND RELAXATION MATRICES
    if (oned)
        cgstiffness;
    else
        cgstiffness2d;
    end
    
    P = floor(P/2);
    PC = floor(P/2);
end

maxlvl = 2;

if (oned)
    kdxloop;
else
    kdxyloop;
end