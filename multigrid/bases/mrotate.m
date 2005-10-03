function [mout] = mrotate(ms)

P = size(ms,1) -1;
temp = ms(2,:);
ms(2:P,:) = ms(3:P+1,:);
ms(P+1,:) = temp;
P = size(ms,2) -1;
temp = ms(:,2);
ms(:,2:P) = ms(:,3:P+1);
ms(:,P+1) = temp;
mout = ms;

return;