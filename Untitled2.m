% A=2;
% D=1/2;
% k=0.4;
% for t=0:1/k:5/k
% 
% x=0:0.1:1;
% u=A*exp(-k*t)*(1-erf(x/(2*sqrt(D*t))));
% hold on
% plot(x,u)
% end
% xlabel('x[L]');
% ylabel('u[M/L]');
% title('Problem 3, (g):');
%%
clc 
clear
A=1;
D=1/2;
c1=2/sqrt(pi);
c2=-2;
x=0:0.1:4;
for t=0:1:5
y=x/sqrt(D*t);
F=c1*(exp(-y.^2/4)+(y/2).*sqrt(pi).*erf(y/2))+(1+c2)*y;
u=sqrt((t*A^2)/D)*F;
hold on
plot(x,u)
end
xlabel('x[L]');
ylabel('u[M/L]');
title('Problem 2, (e):');