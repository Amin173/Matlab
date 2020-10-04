%% Problem 1
%% Linearized Kalman Filter
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

sigma_r=0.05;
sigma_a= 0.4*pi/180;
V=diag([sigma_r,sigma_a,sigma_r,sigma_a].^2);

x_star=[5;0];
P_bar=eye(2);
x_bar=x_star;
Dx_bar=x_bar-x_star;

for i=1:length(T)

z=T(i,2:end)';

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

Dz=z-h;
L=P_bar*H'/(V+H*P_bar*H');

%msmt update
Dx_hat=Dx_bar+L*(Dz-H*(x_bar-x_star));
P_hat=(eye(2)-L*H)*P_bar;
x_hat=x_star+Dx_hat;

%Time update
x_star=[5;1*i]+0.1*randn([2,1]);
P_bar=Phi*P_hat*Phi'+Gamma_w*W*Gamma_w';
Dx_bar=Phi*Dx_hat+Gamma_u*u;
x_bar=x_star+Dx_bar;

    hold on
    plot(i,x_hat(1),'bo')
    plot(i,x_hat(2),'r*')
end

legend('x1Hat [m]','x2Hat [m]')
xlabel('time[s]')
title('Position posterior estimate in a1 and a2 directions (x1 and x2)')
%% Extended Kalman Filter

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

sigma_r=0.05;
sigma_a= 0.4*pi/180;
V=diag([sigma_r,sigma_a,sigma_r,sigma_a].^2);

P_bar=eye(2);
x_bar=[5;0];

for i=1:length(T)

z=T(i,2:end)';

x_star=x_bar;

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

Dz=z-h;
L=P_bar*H'/(V+H*P_bar*H');

%msmt update
Dx_hat=L*Dz;
P_hat=(eye(2)-L*H)*P_bar;
x_hat=x_bar+Dx_hat;

%Time update
P_bar=Phi*P_hat*Phi'+Gamma_w*W*Gamma_w';
x_bar=[5;1*i]+0.1*randn([2,1]);

    hold on
    plot(i,x_hat(1),'bo')
    plot(i,x_hat(2),'r*')
end

legend('x1Hat [m]','x2Hat [m]')
xlabel('time[s]')
title('Position posterior estimate in a1 and a2 directions (x1 and x2)')

%% Iterative Extended Kalman Filter

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

sigma_r=0.05;
sigma_a= 0.4*pi/180;
V=diag([sigma_r,sigma_a,sigma_r,sigma_a].^2);

P_bar=eye(2);
x_bar=[5;0];

eps=1e-3;

for i=1:length(T)

z=T(i,2:end)';


delx=2*eps;
x_star=x_bar;
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

%msmt update
z_star=z-h+H*x_star;
L=P_bar*H'/(V+H*P_bar*H');
x_hat=x_bar+L*(z_star-H*x_bar);
P_hat=(eye(2)-L*H)*P_bar;
delx=((x_hat-x_star)'/P_hat)*(x_hat-x_star)
x_star=x_hat;

end
%Time update

x_bar=[5;1*i]+0.1*randn([2,1]);
P_bar=Phi*P_hat*Phi'+Gamma_w*W*Gamma_w';

    hold on
    plot(i,x_hat(1),'bo')
    plot(i,x_hat(2),'r*')
end
legend('x1Hat [m]','x2Hat [m]')
xlabel('time[s]')
title('Position posterior stimate in a1 and a2 directions (x1 and x2)')


%% Problem 2
clc
clear
syms x1 x2 x3 x4 f w T
F=[0 1 0 0;0 0 (-1/(1-x4)) -(f-x3-w)/(1-x4)^2;0 0 0 0;0 0 0 0];
phi=exp(F*0.1)
%%
clc
clear
Ts=0.1;
Qw=0.05;
Gw=[0;-1;0;0];
Gu=[0;-1;0;0];
H=[1 0 0 0];
V=0.5^2;
P_bar{1}=[0.5^2 0 0 0;0 0.4^2 0 0;0 0 0.1^2 0;0 0 0 (0.2)^2];
gamma_w=Gw*Ts;
W=Qw/Ts;

phi=[ 1, exp(1/10),                    1,                                 1
 1,         1, exp(-0.1), exp(-1)
 1,         1,                    1,                                 1
 1,         1,                    1,                                 1];

gamma_u=[0;0;0;0];

for i=1:10
P=P_bar{i};
L=P*(H'/(V+H*P*H'));
P_hat=(eye(4)-L*H)*P;
P_tupdate=phi*P_hat*phi'+gamma_w*W*gamma_w';
P_bar{i+1}=P_tupdate;

pos_est_err(i,:)=sqrt([P(1,1),P_hat(1,1)]);
scale_est_err(i,:)=sqrt([P(4,4),P_hat(4,4)]);

end

plot(1:10,pos_est_err);
legend('Postion estimate std before measurement: sqrt(Pbar-xx)','Postion estimate std after measurement: sqrt(Phat-xx)');
figure()
plot(1:10,scale_est_err,'r--','bx');
legend('Scale factor estimate std before measurement: sqrt(Pbar-xx)','Scale factor estimate std after measurement: sqrt(Phat-xx)');

%% Problem 4
clc
clear
z=[100;pi/4]+[1*randn([1,1e6]);(pi/180)*randn([1,1e6])];
sigma_r=1;
sigma_th=10*pi/180;
V=diag([sigma_r,sigma_th].^2);

x_bar=[100*cos(pi/4);100*sin(pi/4)];

for i=1:size(z,2)
x1=x_bar(1);
x2=x_bar(2);

d2=x1^2+x2^2;
d=sqrt(d2);
h=[sqrt(x1^2+x2^2);atan2(x2,x1)];
H=[x1/d x2/d;-x2/d2 x1/d2];
P_hat=inv(H'*inv(V)*H);
L=P_hat*H'*inv(V);
x_hat=x_bar+L*(z(:,i)-h);

x1_hat(i)=x_hat(1);
x2_hat(i)=x_hat(2);
end

hold on
plot(x1_hat,x2_hat,'r.')
title('SigmaVOR=10')

plot_gaussian_ellipsoid([100*cos(pi/4);100*sin(pi/4)], P_hat, 2)

legend('Estimated position xHat','Linearized covarience elipse')
xlabel('x1 [km]')
ylabel('x2 [km]')