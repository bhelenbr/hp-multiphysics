elgrid=-1:2./(bsz^(1/(2-oned))):1;
if (oned)
    mode = basissave*evect(bsz*(n-1)+1:bsz*n);
    for ex=0:Nelx-1
        xpos = ex +(elgrid+1.0)/2.0;
        f = double(real(subs(mode*exp(-i*ex*kdx),x,elgrid)));
        plot(xpos,f);
        hold on;
    end
    xlabel('x/\Deltax');
    axis auto;
elseif (Fourier2D)
    %%%%% TO PLOT 2D EIGENFUNCTIONS %%%%
    mode = basissave*evect(bsz*(n-1)+1:bsz*n);
    [egridx,egridy] = meshgrid(elgrid,elgrid);
    for ex=0:Nelx-1
        for ey=0:Nely-1
            f = double(real(subs(mode*exp(-i*(ex*kdx +ey*kdy)),{x,y},{egridx,egridy})));
            xpos = ex +(1.0+egridx)/2.0;
            ypos = ey +(1.0+egridy)/2.0;
            contourf(xpos,ypos,f);
            hold on;
        end
    end
    colorbar
    xlabel('x/\Deltax');
    ylabel('y/\Deltay');
    axis auto;
    axis equal;
else
    %%%%% TO PLOT 2D EIGENFUNCTIONS %%%%
    if (bsz == 1 && maxlvl >1)
        cscale = 2*max(abs(evect));
        hold on
        % SPECIAL CASE FOR AGGLOMERATION
        for ex=0:Nelx-1
            mode = evect(ex*4*bsz*nvar +(n-1)*4*bsz+1:ex*4*bsz*nvar +n*4*bsz);
            for ey=0:Nely-1
                clr = real(mode*exp(-i*ey*kdy))/cscale +0.5;
                
                % WEIRD WAY TO DRAW CONSTANT COLOR LEVEL
                fill([ex*2,ex*2+1,ex*2+1,ex*2],[ey*2,ey*2,ey*2+1,ey*2+1],[clr(1) clr(1) clr(1)]);
                fill([ex*2+1,ex*2+2,ex*2+2,ex*2+1],[ey*2,ey*2,ey*2+1,ey*2+1],[clr(2) clr(2) clr(2)]);
                fill([ex*2+1,ex*2+2,ex*2+2,ex*2+1],[ey*2+1,ey*2+1,ey*2+2,ey*2+2],[clr(3) clr(3) clr(3)]);
                fill([ex*2,ex*2+1,ex*2+1,ex*2],[ey*2+1,ey*2+1,ey*2+2,ey*2+2],[clr(4) clr(4) clr(4)]);
            end
        end
    else    
        [egridx,egridy] = meshgrid(elgrid,elgrid);
        for ex=0:Nelx-1
            mode = basissave*evect(ex*bsz*nvar +(n-1)*bsz+1:ex*bsz*nvar +n*bsz);
            for ey=0:Nely-1
                f = double(real(subs(mode*exp(-i*ey*kdy),{x,y},{egridx,egridy})));
                xpos = ex +(1.0+egridx)/2.0;
                ypos = ey +(1.0+egridy)/2.0;
                contourf(xpos,ypos,f,[-0.45:0.05:0.45]);
                hold on;
            end
        end
        colorbar
        xlabel('x/\Deltax');
        ylabel('y/\Deltay');
        axis auto;
        axis equal;
        myexport(['c' num2str(n)]);
        
        % 1D PLOT FOR SHOWING CONTINUITY
%         figure
%         basis1D = subs(basissave,y,0);
%         for ex=0:Nelx-1
%             mode = basis1D*evect(ex*bsz*nvar +(n-1)*bsz+1:ex*bsz*nvar +n*bsz);
%             g = double(real(subs(mode,x,elgrid)));
%             xpos = ex +(1.0+elgrid)/2.0;
%             plot(xpos,g);
%             hold on;
%         end
%         xlabel('x/\Deltax');
    end

end