function [x,y]=circletheta(r,x0,y0,n,theta0)

x=zeros(1,n);
y=zeros(1,n);


for i =1:n
    theta=pi*(i-1)/(n-1)-theta0;
    [x(i),y(i)]=pol2cart(theta,r);
end

x=x+x0;
y=y+y0;