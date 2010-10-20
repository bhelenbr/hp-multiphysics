% FOR LAGRANGE & SPECTRAL ELEMENT BASIS
if (P > 0)
    [lpt,lwt] = lglnodes(P);
    %lpt = -1:2/P:1
else
  lpt = 0.0;
  lwt = 2.0;
end

if (PC > 0)
  [lptc,lwtc] = lglnodes(PC);
  %lptC = -1:2/PC:1
else
  lptc = 0.0;
  lwtc = 2.0;
end


% THEY COME OUT REVERSE ORDER (URGH!)
[lpt,lsort] = sort(lpt);
lwt = lwt(lsort);
[lptc,lsort] = sort(lptc);
lwtc = lwtc(lsort);

syms x real;


for i1=1:P+1
        basis(i1) = sym('1');
    for i2=1:i1-1
        basis(i1) = simplify(basis(i1)*(x-lpt(i2))/(lpt(i1)-lpt(i2)));
    end
    for i2=i1+1:P+1
        basis(i1) = simplify(basis(i1)*(x-lpt(i2))/(lpt(i1)-lpt(i2)));
    end
end

for i1=1:PC+1
        basisc(i1) = sym('1');
    for i2=1:i1-1
        basisc(i1) = basisc(i1)*(x-lptc(i2))/(lptc(i1)-lptc(i2));
    end
    for i2=i1+1:PC+1
        basisc(i1) = basisc(i1)*(x-lptc(i2))/(lptc(i1)-lptc(i2));
    end
    basisc(i1) = simplify(basisc(i1));
end

% FOR NODAL BASIS
for i1=1:PC+1
    for i2=1:P+1
        restrict(i1,i2) = subs(basisc(i1),x,lpt(i2));
    end
end

matrix;

mapp = diag(lwt);