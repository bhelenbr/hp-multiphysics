
clear
close all

theta = pi/4:0.1:pi
r1 = (cos(theta) + sqrt(cos(theta).^2 + 3*sin(theta).^2/2))./sin(theta).^2
%r2 = (cos(theta) - sqrt(cos(theta).^2 + 3*sin(theta).^2/2))./sin(theta).^2
 
x = r1.*cos(theta)+2;
y = r1.*sin(theta);
plot(x,y)


% syms r theta
% x0 = r*cos(theta)+2
% x1 = r*sin(theta)
% h = x0-1.25-0.5*x1^2
% solve(h,r)

x0 = 10;
x1 = 1.0;
M = 10;
gamma = 1.4;
rhou = 1.4;

uu = 10;
vu = 0;
cu = uu/M;
RTu = cu^2/gamma;
pu = rhou*RTu;


for x1=0.1:4
    r = sqrt((x0-2)^2+x1^2);
    R = 0.5;
    r_exit = sqrt(17.5+64);
    ct_exit = 8/r_exit;

    if ((x0-2)/r < ct_exit) 
        ct = (x0-2)/r ;
    else
        ct = ct_exit;
    end
    st = sqrt(1-ct^2);
    r_wall = sqrt(64+x1^2);
    if (ct > ct_exit) 
        rshock = r_exit;
    else
        if (abs(x1) > 0)
            rshock = (ct + sqrt(ct^2 + 3*st^2/2))/st^2;
        else
            rshock = 0.75;
        end
    end
    norm_mag = 1/((1+(rshock*st)^2)^(1/2));
    norm_x = norm_mag;
    norm_y = -norm_mag*rshock*st;
    vs = 0;
    uu_norm = uu*norm_x+vu*norm_y;
    Md = (uu_norm-vs)/cu;
    pd = pu*(1.0+((2.0*gamma)/(gamma+1.0))*(Md^2-1.0));
    rhod = rhou*(((gamma+1.0)*Md^2)/((gamma-1.0)*Md^2+2.0));
    ud_norm = vs+(rhou/rhod)*(uu_norm-vs);
    RTd = pd/rhod;
    tan_x = -norm_y;
    tan_y = norm_x;
    uu_tan = uu*tan_x+vu*tan_y;
    ud_tan = uu_tan;
    ud = norm_x*ud_norm-norm_y*ud_tan;
    vd = norm_y*ud_norm+norm_x*ud_tan;
    %Velocity varies linearly to zero at cylinder wall = 
    if (ct < ct_exit)
        rfinal = rshock;
    else
        rfinal = r_wall;
    end
    rfinal
    r
    
    ud_use = ud*(r-R)/(rfinal-R)
    vd_use = vd*(r-R)/(rfinal-R);
end


