syms x real

P = input('P');

xdiv = -1:2/(P+1):1

for a1 = 1:P+1
for a2 = 1:P+1
	a(a1,a2) = int(x^(a2-1),xdiv(a1),xdiv(a1+1));
end
pvect(a1) = x^(a1-1);
end


for a1=1:P+1
	rhs(1:P+1) = 0.0;
	rhs(a1) = 2/(P+1);
	bcoeff(a1,:) = (a\rhs')';
	basis(a1) = bcoeff(a1,:)*(pvect');
end

