%% Problem 1
clc
clear
close all
%% (a)
dt=1;
S=1;
phi=[1 dt;0 1];
Q=[S*dt^3/3 S*dt^2/2;S*dt^2/2 S*dt];
H=[1 0];
R=1;
x0=[0;0];
t=1:100;
xt(:,1)=x0;
for i=1:length(t)
    
   w=mvnrnd([0 0],Q);
   z(i)=H*xt(:,i)+normrnd(0,1);
   xt(:,i+1)=phi*xt(:,i)+w';    
   
end
xt(:,length(t)+1)=[];
%% plots
figure(1)
plot(xt(1,:),xt(2,:))
xlabel('x1')
ylabel('x2')
figure(2)
plot(t,z)
xlabel('t[s]')
ylabel('z')
%% (b)
SS=[1;0.01];
alpha=ones([length(SS) 1]);
alpha=alpha/length(SS);
pzs(1)=1;
for i=1:length(SS)
    xhat(:,1,i)=[0;0];
    Phat(:,:,i)=eye(2)*1e4;
end
xa=zeros(2,length(t));
Pa=zeros(2,2,length(t));
for i=1:length(t)
for j=1:length(SS)
    S=SS(j);
    Q=[S*dt^3/3 S*dt^2/2;S*dt^2/2 S*dt];
    
    [x,P]=KFA(xhat(:,i,j),Phat(:,:,j),z(i),Q,R,H,phi);
    
    xhat(:,i+1,j)=x;
    Phat(:,:,j)=P;
    xa(:,i)=alpha(j)*xhat(:,i,j)+xa(:,i);
    Pa(:,:,i)=alpha(j)*P+Pa(:,:,i);
    pz(i+1,j)=pzs(i)*normpdf(z(i),H*x,H*P*H'+R);

end
    alpha=alpha.*pz(i+1,:)';
    alpha=alpha/sum(alpha);
    pzs(i+1)=sum(pz(i+1,:)'.*alpha);
end
for i=1:length(t)
Ct=trace(Pa(:,:,i));
end
%% results
figure(3)
hold on
plot(xt(1,:),xt(2,:));
plot(xa(1,:),xa(2,:),'--');
legend('xt','xa')
xlabel('x1')
ylabel('x2')
hold off

figure(4)
hold on
plot(t,z)
plot(t,H*xt,'o-')
xlabel('t[s]')
legend('Observation','Filtered Observation')

figure(5)
plot(t,Ct,'*')
xlabel('t [s]')
ylabel('AKF Cov. trace')
%% Problem 2
clc 
clear
close all
%% Importing data
T=importdata('C:\Users\amink\Documents\MATLAB\T_sensor_data.mat');
H=importdata('H.mat');
%% Parameters Initialization
m=100;
Sw=0.1; %process noise
Sv=0.1; %measurement noise
V=eye(size(T,2))*Sv^2;
W=eye(m)*Sw^2;
gw=ones([m 1])';
phi=zeros([m,m]);
phi(1,2)=1;
phi(m,m-1)=1;
for i=2:m-1
    phi(i,i+1)=1/2;
    phi(i,i-1)=1/2;
end

H=zeros([size(T,2) m]);
mz=10:10:90;
for i=1:length(mz)
    H(i,mz(i))=1;
end
%% Kalman Filter
x_bar=0.01*((1:m)-50).^2;
x_bar=reshape(x_bar,[m,1]);
P_bar=eye(m)*Sv^2;
for i=1:size(T,1)
    L=P_bar*H'*inv(H*P_bar*H'+V);
    x_hat=x_bar+L*(T(i,:)'-H*x_bar);
    P_hat=P_bar-L*H*P_bar;
    P_bar=phi*P_hat*phi'+gw*W*gw';
    x_bar=phi*x_hat;
    
    x_5_sav(i)=x_hat(5);
    x_45_sav(i)=x_hat(45);
    P_5_sav(i)=P_hat(5,5);
    P_45_sav(i)=P_hat(45,45);
    
end
%% result
dt=1;
t=1:dt:size(T,1);

figure(6)
plot(t,x_5_sav,':',t,x_45_sav,'-')
xlabel('time [s]')
legend('That(j=5)','That(j=45)')
ylabel('Temp [^oC]')
title('Kalman Filter Solution')

figure(7)
plot(t,P_5_sav,':',t,P_45_sav,'-')
xlabel('time [s]')
legend('Phat(5,5)','Phat(45,45)')
ylabel('Temp [^oC]')
title('Kalman Filter Solution - Temp. Covariance')
%% Emsemble Kalman Filter
clc
clear
close all
%% Importing data
T=importdata('C:\Users\amink\Documents\MATLAB\T_sensor_data.mat');
m=100;
Sw=0.1; %process noise
Sv=0.1; %measurement noise
V=eye(size(T,2))*Sv^2;
W=eye(m)*Sw^2;
gw=ones([m 1])';
phi=zeros([m,m]);
phi(1,2)=1;
phi(m,m-1)=1;
for i=2:m-1
    phi(i,i+1)=1/2;
    phi(i,i-1)=1/2;
end

H=zeros([size(T,2) m]);
mz=10:10:90;
for i=1:length(mz)
    H(i,mz(i))=1;
end

x_bar=0.01*((1:m)-50).^2;
x_bar=reshape(x_bar,[m,1]);
P_bar=eye(m)*Sv^2;
%%
dN=100;
N=floor(1:100/dN:100);  %index of sampled states
M=100;        %number of sampling points for each state
Xbar_i=zeros(M,length(N));
for i=1:length(N)
Xbar_i(:,i)=normrnd(x_bar(N(i)),P_bar(N(i),N(i)),[M 1]);
end
% for i=1:size(T,1)
% Z(:,i)=T(i,:)'+normrnd(0,Sv,[size(T,2) 1]); 
% end
for i=1:size(T,1)
Z=T(i,:)'+normrnd(0,Sv,[size(T,2) length(N)]); 
Xbar=mean(Xbar_i,2);
EXbar=Xbar_i-Xbar;
P_bar=EXbar*EXbar'/(length(N)-1);
S=(H*EXbar)*(H*EXbar)'/(length(N)-1);
L=EXbar*(H*EXbar)'*inv(S+V)/(length(N)-1);
Xhat=Xbar_i+L*(Z-H*Xbar_i);
Phat=EXbar*EXbar'/(length(N)-1);
W=normrnd(0,Sw,[m length(N)]);
Xbar_i=phi*Xhat+gw'.*W;
xhat=mean(Xhat,2);
    x_5_sav(i)=xhat(5);
    x_45_sav(i)=xhat(45);
    P_5_sav(i)=Phat(5,5);
    P_45_sav(i)=Phat(45,45);
end
%% result
dt=1;
t=1:dt:size(T,1);

figure(6)
plot(t,x_5_sav,':',t,x_45_sav,'-')
xlabel('time [s]')
legend('That(j=5)','That(j=45)')
ylabel('Temp [^oC]')
title('Ensemble Kalman Filter Solution')

figure(7)
plot(t,P_5_sav,':',t,P_45_sav,'-')
xlabel('time [s]')
legend('Phat(5,5)','Phat(45,45)')
ylabel('Temp [^oC]')
title('Ensemble Kalman Filter Solution - Temp. Covariance')