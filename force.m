clc
clear
F1=100;
F2=10;
theta=pi/4;
x=1:0.1:10;
f=1000;
w=2*pi*f;
fun1 = @(x)F1*[cos(theta)].*cos(w*x)-F2*[cos(theta)].*sin(w*x);
fun2 = @(x)F1*[sin(theta)].*cos(w*x)-F2*[sin(theta)].*sin(w*x);
vx=integral(fun1,0,10);
vy=integral(fun2,0,10);
x=integral(vx,0,10);
y=integral(vy,0,10);

% subplot(3,1,1)
% plot(fit1,t,) % cfit plot method
% subplot(3,1,2)
% plot(xdata,d1,'m') % double plot method
% grid on
% legend('1st derivative')
% subplot(3,1,3)
% plot(xdata,d2,'c') % double plot method
% grid on
% legend('2nd derivative')
%% 
clc
clear



syms t F1 F2 theta f
w=2*pi*f;
fun1 = F1*[cos(theta) sin(theta)].*cos(w*t)-F2*[cos(theta) sin(theta)].*sin(w*t);
v=int(fun1,t);
x=simplify(int(v,t))
%%
clc
clear
t=[0:0.1:100]';
theta=pi/4;
f=1000;
F1=0;
F2=100;
x12=[ -(cos(theta)*(F1 + F1*cos(2*pi*f*t) - F2*sin(2*pi*f*t) - 2*pi*F2*f*t))/(4*f^2*pi^2), -(sin(theta)*(F1 + F1*cos(2*pi*f*t) - F2*sin(2*pi*f*t) - 2*pi*F2*f*t))/(4*f^2*pi^2)];

plot(t,x12(:,1));
