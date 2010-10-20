
function [krtn]=create_operator(K,lvl,nb)
    global oned Fourier2D kdx kdy Nelx
    mid = floor(nb/2) +1;
    krtn = zeros(size(K{lvl,mid,mid}));
    if (oned)
        for m = 1:nb
            krtn = krtn +K{lvl,mid,m}*exp(-i*(m-mid)*kdx);
        end
    elseif (Fourier2D)
        for m=1:nb
            for n=1:nb
                krtn = krtn +K{lvl,m,n}*exp(-i*((m-mid)*kdy +(n-mid)*kdx));
            end
        end
    else
        % FOURIER IN Y ONLY
        krtn = blktimes(zeros(Nelx,Nelx),K{lvl,mid,mid});
        for n=1:nb
            ktemp = zeros(size(K{lvl,mid,mid}));
            for m=1:nb
                ktemp = ktemp+ K{lvl,m,n}*exp(-i*((m-mid)*kdy));
            end

            pos = n-mid;
            krtn = krtn + blktimes(diag(ones(Nelx-abs(pos),1),pos),ktemp);
        end
    end
end