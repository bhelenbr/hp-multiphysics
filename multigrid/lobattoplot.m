N = 20;
for i = 2:N
[wgt,x] = lobatto(i);
dx(i-1) = abs(x(1)-x(2));
p(i-1) = i;
end

plot(p,dx);
figure
loglog(p,dx,'x');
xlabel('Number of Gauss-Lobatto Points');
ylabel('Minimum Point Spacing');