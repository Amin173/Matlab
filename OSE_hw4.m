%% Problem 1

clc
clear
x0=0;
P0=1;
x_bar(1)=x0;
P_bar(1)=P0;
Phi=1;
V=1;
W=0;
Gamma_w=0;
Gamma_v=1;
z=[-1.43 -2.6656 -0.7123 -2.1466 0.1909 0.1892 -1.0376 -0.6727 -0.8254];
for i=1:length(z)
    L(i)=P_bar(i)*Phi*inv(V+Phi^2*P_bar(i));
    x_hat(i)=x_bar(i)+L(i)*(z(i)-Phi*x_bar(i));
    P_hat(i)=(1-L(i)*Phi)*P_bar(i);
    P_bar(i+1)=Phi^2*P_hat(i);
    x_bar(i+1)=Phi*x_hat(i);
end
figure(1)

plot(1:length(z),x_hat,1:length(z),mean(z)*z./z);
xlabel('Measurment No.')
legend('xhat','zbar')
figure(2)
sigma=P_hat.^0.5;
plot(1:length(z),sigma)
xlabel('Measurment No.')
ylabel('sigma_x')

%% Problem 2

clc
clear
V=0.15^2;
W=0.03*V;
N=20;
alpha=0.9;
x0=1;
sigma0=0.1^2;

w=rand([N,1])*W^0.5;
v=rand([N,1])*V^0.5;

x(1)=x0;
P_bar(1)=sigma0;
x_bar(1)=x0;

for i=1:N
    z(i)=x(i)+v(i);
    x(i+1)=alpha*x(i)+w(i);
end

for i=1:N
    L(i)=P_bar(i)/(V+P_bar(i));
    x_hat(i)=x_bar(i)+L(i)*(z(i)-x_bar(i));
    P_hat(i)=(1-L(i))*P_bar(i);
    P_bar(i+1)=alpha^2*P_hat(i)+W;
    x_bar(i+1)=alpha*x_hat(i);
end

figure(1)
plot(1:N,x_hat,1:N,x_hat-x(1:N))
legend('State Estimate(xhat)','State Estimate Error')
xlabel('Measurement No.')
title('Post measurement state estimation')
figure(2)
plot(1:N,x_bar(1:N),1:N,x_bar(1:N)-x(1:N))
legend('State Estimate(xbar)','State Estimate Error')
xlabel('Measurement No.')
title('Pre-measurement state estimation')
figure(3)
plot(1:N,sqrt(P_hat),1:N,sqrt(P_bar(1:N)))
legend('sqrt(Phat)','sqrt(Pbar)')
xlabel('Measurement No.')

%% Problem 3

clc
clear
Q=0.01;
T=0.1;
V=0.2^2;
[phi,gamma_u]=c2d([0 1;0 0],[0;0],T);
H=[1 0];
Gamma_w=[0;0.1];
W=Q/T;
N=11;


P_bar=1e9*eye(2);
for i=1:N
    P_hat=inv(inv(P_bar)+H'/V*H);
    P_bar=phi*P_hat*phi'+Gamma_w*W*Gamma_w';
    Px(i)=P_bar(1,1);
    Px_hat(i)=P_hat(1,1);
end
hold on
plot(sqrt(Px(2:end)))
plot(sqrt(Px_hat(2:end)))
xlabel('Time[T]')
legend('Pre-measurement variance','Post-measurement variance')
%% Problem 4
clc 
clear
syms x1 x2
xv=[x1 x2];
xa=[0 0];
xb=[10 10];
h=[norm(xv-xa);atan(x1/x2);norm(xv-xb);atan((10-x1)/(10-x2))];
[simplify(diff(h,x1)) simplify(diff(h,x2))]

%%
clc
clear
T = importdata('las_hw.dat');
T(:,[3 5])=T(:,[3 5])*pi/180;
T(:,3)=T(:,3)-pi;
F=zeros(2,2);
Gu=[0;1];
u=1;
[Phi,Gamma_u]=c2d(F,Gu,1);
Gw=[1 1];
W=0.1^2;
Gamma_w=Gw*1;

P_bar=eye(2);
x_bar=[5;0];
sigma_r=0.05;
sigma_a= 0.4*pi/180;
V=diag([sigma_r,sigma_a,sigma_r,sigma_a].^2);
for i=1:length(T)
    hold on
    plot(i,x_bar(1),'bo')
    plot(i,x_bar(2),'r*')
z=T(i,2:end)';
x_star=x_bar;
eps=1e-3;
delx=2*eps;
while delx>eps
x1=x_star(1);
x2=x_star(2);
h=[((x1)^2 + (x2)^2)^(1/2);
   atan2(x1,x2);
 ((x1 - 10)^2 + (x2 - 10)^2)^(1/2);
  atan2(-(x1 - 10),-(x2 - 10))];
           
H=[x1/(x1^2 + x2^2)^(1/2),  x2/(x1^2 + x2^2)^(1/2)
   x2/(x1^2 + x2^2),  -x1/(x1^2 + x2^2)
  (x1-10)/((x1 - 10)^2 + (x2 - 10)^2)^(1/2),  (x2-10)/((x1 - 10)^2 + (x2 - 10)^2)^(1/2)
   1/(((x1 - 10)^2/(x2 - 10)^2 + 1)*(x2 - 10)),   -(x1 - 10)/(((x1 - 10)^2/(x2 - 10)^2 + 1)*(x2 - 10)^2)];

z_star=z-h+H*x_star;
L=P_bar*H'/(V+H*P_bar*H');
x_hat=x_bar+L*(z_star-H*x_bar);
P_hat=(eye(2)-L*H)*P_bar;
delx=((x_hat-x_star)'/P_hat)*(x_hat-x_star);
x_star=x_hat;
end
P_bar=Phi*P_hat*Phi'+Gamma_w*W*Gamma_w';
x_bar=Phi*x_hat+Gamma_u*u;
end
legend('x1bar [m]','x2bar [m]')
xlabel('time[s]')
title('Position Estimate in a1 and a2 directions (x1 and x2)')