function [correction] = mgcycle2d(vw,lvl,rsdlin)
global KB MB RB vrestrict vmaxeig;
global kdx kdy omega maxlvl sys_flag rlx_flag sweep_flag;

    % LOAD MATRICES
	bsz = size(KB{lvl,3,3},1);
    
	% STIFFNESS
	k1 = zeros(bsz,bsz);
	for m=1:5
		for n=1:5
            k1 = k1 +KB{lvl,m,n}*exp(i*((m-3)*kdy +(n-3)*kdx));
		end
	end
        
    if (lvl == maxlvl)
        % SOLVE MATRIX DIRECTLY 
    	correction = inv(k1)*rsdlin;
    else
        % PRECONDITIONER
        % maxeig = vmaxeig(lvl);
        
	    % PRECONDITIONER
        if (sweep_flag ~= 2 && sweep_flag ~= 5)
            m1 = zeros(bsz,bsz);
            for m=1:5
                for n=1:5
                    m1 = m1 +MB{lvl,m,n}*exp(i*((m-3)*kdy +(n-3)*kdx));
                end
            end
            m1 = inv(m1);
        elseif (sweep_flag == 2)
            mbl = zeros(bsz,bsz);
            for m=1:3
                for n=1:3
                    mbl = mbl +MB{lvl,m,n}*exp(i*((m-3)*kdy +(n-3)*kdx));
                end
            end
            mtr = zeros(bsz,bsz);
            for m=3:5
                for n=3:5
                    mtr = mtr +MB{lvl,m,n}*exp(i*((m-3)*kdy +(n-3)*kdx));
                end
            end   
            m1 = inv(mtr) +inv(mbl) -inv(mtr)*k1*inv(mbl);
        else
            mb = zeros(bsz,bsz);
            for m=1:3
                for n=1:5
                    mb = mb +MB{lvl,m,n}*exp(i*((m-3)*kdy +(n-3)*kdx));
                end
            end
            mt = zeros(bsz,bsz);
            for m=3:5
                for n=1:5
                    mt = mt +MB{lvl,m,n}*exp(i*((m-3)*kdy +(n-3)*kdx));
                end
            end   
            m1 = inv(mt) +inv(mb) -inv(mt)*k1*inv(mb);
        end
        
        % RESTRICTION OPERATOR
        bszc = size(RB{lvl,3,3},1);
        r1 = zeros(bszc,bsz);
        for m=1:3
            for n=1:3
                r1 = r1 +RB{lvl,m,n}*exp(i*((m-2)*kdy +(n-2)*kdx));
            end
        end
        
        if (lvl == 1)
            rsdlin = zeros(bsz,bsz);
            correction = eye(bsz,bsz);
        else
            correction = zeros(size(rsdlin));
        end
        
        for cyc = 1:vw
            
            % FINE GRID AMPLIFICATION FACTOR
            if (sys_flag == 1 || rlx_flag == 2)
                correction = correction -m1*(k1*correction -rsdlin);
            else
                a = m1*k1;
		        correction = rk3_5(correction,zeros(size(a)),2.0*a,-m1*rsdlin);
            end
            
            cchange = mgcycle2d(vw,lvl+1,r1*(k1*correction -rsdlin));
            
            % CONJUGATE TRANSPOSE OF RESTRICTION FOR CONTINUOUS
            % DOESN'T MATTER FOR DISCONTINUOUS
            correction = correction -r1'*cchange;
        end
    end
return;