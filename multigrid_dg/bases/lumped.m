clear restrict
clear ifmv ifmvc CM RW a b c;

dubiner;


% FOR LUMPING BASIS
% FIND FINE GRID LUMPED MASS MATRIX
ifmv = zeros(P+1,P+1);
if (P > 1) 
	CM = [2:P-1,P+1];
	RW = [2:P];
	b = ms(1,CM)';
	a = ms(RW,CM)';
	c = a\b;
	ifmv(1,RW) = c';
	CM = [1,2:P-1];
	RW = [2:P];
	b = ms(P+1,CM)';
	a = ms(RW,CM)';
	c = a\b;
	ifmv(P+1,RW) = c';
    ifmv = eye(P+1,P+1) -ifmv;
else
	ifmv = eye(P+1,P+1);
end

ifmvc = zeros(PC+1,PC+1);
if (PC > 1) 
	CM = [2:PC-1,PC+1];
	RW = [2:PC];
	b = msc(1,CM)';
	a = msc(RW,CM)';
	c = a\b;
	ifmvc(1,RW) = c';
	CM = [1,2:PC-1];
	RW = [2:PC];
	b = msc(PC+1,CM)';
	a = msc(RW,CM)';
	c = a\b;
	ifmvc(PC+1,RW) = c';
    ifmvc = eye(PC+1,PC+1) -ifmvc;
else
	ifmvc = eye(PC+1,PC+1);
end

basis = (ifmv*(basis'))';


restrict = ifmvc*restrict*inv(ifmv);


matrix;

% ms = ifmv*ms;
% cv = ifmv*cv;
% df = ifmv*df;

% mapp = ms;
% if (P > 1)
%     %mapp(P,P) -= mapp(P,1)/mapp(1,1)*mapp(1,P);
% 	mapp(1,P) = 0.0;
% 	%mapp(P,P) -= mapp(P,P+1)/mapp(P+1,P+1)*mapp(P+1,P);
% 	mapp(P+1,P) = 0.0;
% else
%     mapp(1,2) = 0.0;
% 	mapp(2,1) = 0.0;
% end
% eig(inv(mapp)*(mapp-ms))



