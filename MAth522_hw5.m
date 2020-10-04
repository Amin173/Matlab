clc;clear
A=3;
k=0.693/5;
t=1:0.01:8;
x=(A/k)*((1-exp(-k*t))-(1-exp(-k*(t-2))).*((t>2)+0.5*(t==2)));
plot(t,x);
ylabel('x[M]')
xlabel('time[hr]')
title('Problem 1(e)')
%%
clc 
clear
A=400;
k=0.9;
k2=0.693/2.3;
t=1:0.2:8;
x=(A)*((exp(-k*t))+(exp(-k*(t-1/2))).*((t>1/2)+0.5*(t==1/2))+(exp(-k*(t-3/2))).*((t>3/2)+0.5*(t==3/2))+(exp(-k*(t-5/2))).*((t>5/2)+0.5*(t==5/2))+(exp(-k*(t-7/2))).*((t>7/2)+0.5*(t==7/2))+(exp(-k*(t-9/2))).*((t>9/2)+0.5*(t==9/2)));
y=(k2*A/(k2-k))*((exp(-k*t)-exp(-k2*t))+(exp(-k*(t-1/2))-exp(-k2*(t-1/2))).*((t>1/2)+0.5*(t==1/2))+(exp(-k*(t-3/2))-exp(-k2*(t-3/2))).*((t>3/2)+0.5*(t==3/2))+(exp(-k*(t-5/2))-exp(-k2*(t-5/2))).*((t>5/2)+0.5*(t==5/2))+(exp(-k*(t-7/2))-exp(-k2*(t-7/2))).*((t>7/2)+0.5*(t==7/2))+(exp(-k*(t-9/2))-exp(-k2*(t-9/2))).*((t>9/2)+0.5*(t==9/2)));
plot(t,x,'-s',t,y,'-o');
legend('Compartment 1','Compartment 2')s
ylabel('x[mg]')
xlabel('time[day]')
title('Problem 2(b)')