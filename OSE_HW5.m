
%% Problem 1
clc
clear

Ts=0.1;
Qw=0.05;
Gw=[0;-1;0];
Gu=[0;1;0];
F=[0 1 0;0 0 -1;0 0 0];
H=[1 0 0];
V=0.5^2;
P_bar{1}=[0.5^2 0 0;0 0.4^2 0;0 0 0.1^2];
[phi,gamma_u]=c2d(F,Gu,Ts);
gamma_w=Gw*Ts;
W=Qw/Ts;

for i=1:10
P=P_bar{i};
L=P*(H'/(V+H*P*H'));
P_hat=(eye(3)-L*H)*P;
P_tupdate=phi*P_hat*phi'+gamma_w*W*gamma_w';
P_bar{i+1}=P_tupdate;

pos_est_err(i,:)=sqrt([P(1,1),P_hat(1,1)]);
end

plot(1:10,pos_est_err);
legend('Postion estimate std before measurement: sqrt(Pbar-xx)','Postion estimate std after measurement: sqrt(Phat-xx)');
xlabel('Time step (Ts=0.1[s])');

%% Problem 3
clc
clear
F=[0 1;0 0];
Gu=[0;0];
Gw=[0;1];
H=[1 0];
V=4;
Ts=0.1;
Qw=0.01;
[phi,gamma_u]=c2d(F,Gu,Ts);
gamma_w=Gw*Ts;
W=Qw/Ts;

P_b=eye(2)*1e9;

for i=1:100
    P_b=phi*P_b*phi'+gamma_w*W*gamma_w';
    P_bar{i+1}=P_b;
    if rem(i,10)==0
        P_h=inv(inv(P_b)+H'/V*H);
        P_bar{i+1}=P_h;
        P_b=P_h;
    end
    sigma_h(i)=P_b(1,1)^0.5;
end
hold on
plot(1:100,sigma_h);
xlabel('Time steps (0.1[s])');
ylabel('Position error standard dev.');
title('Measurements every 1[s]');

for i=101:200
    P_b=phi*P_b*phi'+gamma_w*W*gamma_w';
    P_bar{i+1}=P_b;
%     if rem(i,10)==0
%         P_h=inv(inv(P_b)+H'/V*H);
%         P_bar{i+1}=P_h;
%         P_b=P_h;
%     end
    sigma_h(i)=P_b(1,1)^0.5;
end
figure()
plot(50:200,sigma_h(50:200));
xlabel('Time steps (0.1[s])');
ylabel('Position error standard dev.');
title('No measurements after t=10[s]');

%% Problem 4
clc
clear

%Parameters
alpha=0.95;
phi=alpha;
gamma_w=1;
H=1;
gamma_v=1;


Sv=0.2;
Sw=0.07;
W=Sw^2;
V=Sv^2;
N=20;      %No. of time sampling
w=Sw*randn(N,1);
v=Sv*randn(N,1);



%Initial Estimate
x_bar(1)=1;
P_bar(1)=0.6^2;

%Real position
x(1)=x_bar(1)+sqrt(P_bar)*randn(1);

for i=1:N
    
    % Real position propagation
    x(i+1)=phi*x(i)+w(i);
    
    % measurement generation 
    z(i)=H*x(i)+v(i);
    % Kalman gain
    L=P_bar*H/(V+H*P_bar*H');
    % Time update
    P_hat(i)=(1-L*H)*P_bar(i);
    x_hat(i)=x_bar(i)+L*(z(i)-H*x_bar(i));
    
    %Time update
    x_bar(i+1)=phi*x_hat(i);
    P_bar(i+1)=phi*P_hat(i)*phi'+gamma_w*W*gamma_w';

    
end


% Batch smoothing process
zB(1)=x_bar(1);
vB(1)=P_bar(1)*randn(1);
VB(1)=P_bar(1);
for i=1:N
    zB(2*i)=z(i);
    zB(2*i+1)=0;
    vB(2*i)=V*randn(1);
    vB(2*i+1)=W*randn(1);
    VB(2*i,2*i)=V;
    VB(2*i+1,2*i+1)=W;
end

HB=zeros(2*N+1,N);
HB(1:3,1)=[1;1;phi];
   
    for i=2:N
        HB(:,i)=[0;0;0;HB(1:(2*N+1)-3,i-1)];
    end
    
PB_bar=diag(P_hat);
xB_bar=x_hat';

L=PB_bar*HB'/(VB+HB*PB_bar*HB');
xB_hat=xB_bar+L*(zB'-HB*xB_bar);
PB_hat=(eye(N)-L*HB)*PB_bar;


t=1:20;
hold on
plot(t,x(t));
errorbar(t,x_bar(t),sqrt(P_bar(t)));
errorbar(t,x_hat(t),sqrt(P_hat(t)));
errorbar(t,xB_hat,sqrt(diag(PB_hat)));
xlabel('Time [s]');
legend('Real position x(t)','Prior estimate xbar','Post measurement estimate xhat','Smoothed data');


