%% hw3 problem 2
clc; clear;
%solution using matlab
k=[0.5,1,2];
for i=1: 3
num=k(i);
den=[1 k(i)-1];
T=1;
Y=tf(num,den,T); 
subplot(2,2,i)
hold on
step(Y,10)
t=0:T:10;
yt=(k(i)-1).^t+(t>1);
title(['K=',num2str(k(i))])
plot(t,yt,'--')
legend('Matlab solution','Solution by hand')
end

%% hw3 problem3
clc
clear
T=1/10;
num=[1 0];
den=[1 -2 1];
sys1=tf(num,den,T);
sys2=tf([1 0],[1 -1],T);
t=0:T:10*T;
% u=t/T;
% u=sin(0.1*ws*t/T);
u=ones(length(t),1);
% step(sys)

subplot(2,2,1)
hold on
lsim(sys1,t,t,'-sk',2);
lsim(sys2,t,t,'-og',2);
legend('FOH','ZOH')
title('ramp function input')

subplot(2,2,2)
hold on
lsim(sys1,u,t,'-sk',2);
lsim(sys2,u,t,'-og',2);
legend('FOH','ZOH')
title('unit step function input')

subplot(2,2,3)
hold on
lsim(sys1,sin((2*pi/T)*t),t,'-sk',2)
lsim(sys2,sin((2*pi/T)*t),t,'-og',2)
legend('FOH','ZOH')
title('sin function input')


%% hw3 problem4
clc
clear

T=0.25;
% sys=tf([1 1],[0.1 1]);
% opt = c2dOptions('Method','foh');
% dsys=c2d(sys,T,opt)
dsys=tf([1+T -1],[T+0.1 -0.1],T)
[mag,phase]=bode(dsys,3)


