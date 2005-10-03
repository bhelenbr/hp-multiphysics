function [correction] = mgcycle(vw,lvl,rsdlin)
global KB MB RB vrestrict;
global omega kdx maxlvl sys_flag rlx_flag sweep_flag nvar;

	% STIFFNESS
    k1 = zeros(size(KB{lvl,3}));
    for m = 1:5
	    k1 = k1 +KB{lvl,m}*exp(-i*(m-3)*kdx);
    end
    
    if (lvl == maxlvl)
        % SOLVE MATRIX DIRECTLY 
    	correction = inv(k1)*rsdlin;
    else
		restrict = vrestrict{lvl};
		bsz = size(KB{lvl,3},1);
        
        % RELAXATION
        if (sweep_flag ~= 2) 
			m1 = zeros(size(KB{lvl,3}));
			for m = 1:5
                m1 = m1 + MB{lvl,m}*exp(-i*(m-3)*kdx);
			end
            m1 = inv(m1);
        else
            ml = zeros(size(KB{lvl,3}));
			for m = 1:3
                ml = ml + MB{lvl,m}*exp(-i*(m-3)*kdx);
			end   
            mr = zeros(size(KB{lvl,3}));
			for m = 3:5
                mr = mr + MB{lvl,m}*exp(-i*(m-3)*kdx);
			end   
            m1 = inv(mr) +inv(ml) -inv(mr)*k1*inv(ml);
        end
        
        % RESTRICTION OPERATOR
		r1 = zeros(size(RB{lvl,2}));
		for m = 1:3
            r1 = r1 + RB{lvl,m}*exp(-i*(m-2)*kdx);
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
		        correction = rk3_5(correction,zeros(size(a)),+2.0*a,-m1*rsdlin);
            end

            cchange = mgcycle(vw,lvl+1,r1*(k1*correction -rsdlin));
            
            % CONJUGATE TRANSPOSE OF RESTRICTION FOR CONTINUOUS
            % DOESN'T MATTER FOR DISCONTINUOUS
            correction = correction -r1'*cchange;
        end
    end
return;