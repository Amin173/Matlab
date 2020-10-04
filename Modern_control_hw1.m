clc; clear;
%solution using matlab
num=[2 -2];
den=[1 -3 2];
T=0.1;
t=0:T:100*T;
u=t/T;

Y=tf(num,den,0.1); 
lsim(Y,u,t);

%solution by hand
hold on
k=0:100;
y=(2.^(k+1));
plot(k*T,y,':','linewidth',2);

legend('Matlab solution','Direct solution (by hand)')
% 
%%
% clc
% clear
% T=0.5;
% sys=tf(3,conv([1 1],[1 3]));
% H=c2d(sys,T,'foh');
% pzmap(H)
% grid on
