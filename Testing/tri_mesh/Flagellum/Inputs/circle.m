function [x,y]=circle(r,x0,y0,n)

x=zeros(1,n);
y=zeros(1,n);


for i =1:n
    theta=2*pi*(i-1)/(n-1);
    [x(i),y(i)]=pol2cart(theta,r);
end

x=x+x0;
y=y+y0;