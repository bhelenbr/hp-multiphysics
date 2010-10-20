function [correction,rka] = mgcycle2d(vw,lvl,rsdlin)
global KB MB MB1 RB D Ly Uy Lx Ux omega;
global kdx kdy maxlvl analysis_type rk_flag oned Fourier2D Nelx rlx_flag

	% STIFFNESS
    k1 = create_operator(KB,lvl,5);
        
    if (lvl == maxlvl && maxlvl > 1)
        % SOLVE MATRIX DIRECTLY 
    	correction = k1\rsdlin;
    else
        % PRECONDITIONER
        if (analysis_type == 0)
            % ONE SWEEP SCHEMES
            m1 = create_operator(MB,lvl,5);
            a = m1\k1;    
        elseif (analysis_type == 1)
            % TWO SWEEP SCHEMES
            mbl = create_operator(MB,lvl,5);
            mtr = create_operator(MB1,lvl,5);
            
            % TWO STEP INVERSION WITH RESIDUAL RE-EVALUATION %
            a= -((eye(size(mtr))-mtr\k1)*(eye(size(mbl))-mbl\k1)-eye(size(mtr)));
        
        elseif (analysis_type == 2)
            
            if (rlx_flag ==7 && Fourier2D == 0)
                %ILU
                %LOWER AND UPPER TRIANGULAR PARTS
                bsz=size(KB{lvl,3,3},1);
                md = eye(bsz,bsz)/abs(omega); % other components of L are scaled by omega in dgstiffness2d
                for m = 1:Nelx                     
                    %LOWER TRIANGULAR DIAGONAL
                    mdiag{m} = md + Ly{lvl,m}*exp(-i*kdy);

                    %UPPER TRIANGULAR DIAGONAL
                    mdiag1{m} = D{lvl,m} + Uy*exp(i*kdy);
                end
                    
                % FULL LOWER TRIANGULAR MATRIX
                % Put diagonal blocks in
                mbl = blkdiag(mdiag{:});
                % Make subdiagonal)
                subdiag = blkdiag(Lx{lvl,:});
                % don't forget to clear all otherwise weird shit happens
                % Put bsz zeroz on top
                subdiag = vertcat(zeros(bsz,bsz*(Nelx-1)), subdiag);
                % Put bsz zeros on right
                subdiag = horzcat(subdiag,zeros(bsz*Nelx,bsz));
                mbl = mbl + subdiag;

                % FULL UPPER TRIANGULAR MATRIX
                mtr = blkdiag(mdiag1{:}) + blktimes(diag(ones(Nelx-1,1),1),Ux);
                
            else
                
                mbl = create_operator(MB,lvl,5);
                mtr = create_operator(MB1,lvl,5);
                
            end
           
            % TWO STEP INVERSION WITH NO RESIDUAL REEVALUATION 
            a = mtr\(mbl\k1);
        end       
        
        % RESTRICTION OPERATOR
        r1 = create_operator(RB,lvl,3);
        
        if (lvl == 1)
            rsdlin = zeros(size(k1));
            correction = eye(size(k1));

            % FINE GRID AMPLIFICATION FACTOR
            if (~rk_flag)
                rka = eye(size(a)) -a;
            else
                rka = rk3_5(eye(size(a)),zeros(size(a)),a,zeros(size(a)));
            end
        else
            correction = zeros(size(rsdlin));
        end
        
        % IF maxlvl = 1 vw gets set to 0 to skip this
        for cyc = 1:vw
            
            % FINE GRID AMPLIFICATION FACTOR
            if (~rk_flag)
                correction = correction -a*(correction -k1\rsdlin);
            else
                correction = rk3_5(correction,zeros(size(a)),a,-a/k1*rsdlin);
                % correction = rk3_5(correction,zeros(size(a)),a,-a/k1*rsdlin);
            end
            
            cchange = mgcycle2d(vw,lvl+1,r1*(k1*correction -rsdlin));
            
            % CONJUGATE TRANSPOSE OF RESTRICTION FOR CONTINUOUS
            % DOESN'T MATTER FOR DISCONTINUOUS
            correction = correction -r1'*cchange;
            
            % FINE GRID AMPLIFICATION FACTOR
            % if (~rk_flag)
            %     correction = correction -a*(correction -k1\rsdlin);
            % else
            %     correction = rk3_5(correction,zeros(size(a)),a,-a/k1*rsdlin);
            %     correction = rk3_5(correction,zeros(size(a)),a,-a/k1*rsdlin);
            % end
        end
    end
end
    
    
    
        
            
