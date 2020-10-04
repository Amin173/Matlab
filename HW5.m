clear
clc
y0=[-45*pi/180;45*pi/180;0;0];
tf=10;
global kp kd I1 I2 
kp=600;
kd=500;
[t,x]=ode45(@equ,[0 tf],y0);

% 
% A=xy(x(:,1),x(:,2));
% hold on
% plot(A(:,1),A(:,2));
% xlabel('X [m]')
% ylabel('Y [m]')
% title('End effector trajectory')
% txt1 = '\leftarrow start point';
% txt2='\leftarrow end point';
% text(A(1,1),A(1,2),txt1)
% text(A(end,1),A(end,2),txt2)
% hold on 
a=2*sin(2*t);
b=sin(t);
% B=xy(a,b);
% hold on
% plot(B(:,1),B(:,2),'k:','linewidth',2);
% legend('Followed Path','Desired path')

% subplot(2,2,1)
% plot(t,x(:,1));
% xlabel('Time [s]')
% ylabel('Joint Angle [rad]')
% title('Joint 1 angular position')
% hold on
% plot(t,a,'k:','linewidth',2)
% legend('PD controler output','Desired angle')
% 
% subplot(2,2,2)
% plot(t,x(:,2));
% xlabel('Time [s]')
% ylabel('Joint Angle [rad]')
% title('Joint 2 angular position')
% hold on
% plot(t,b,'k:','linewidth',2)
% legend('PD controler output','Desired angle')
% 
% subplot(2,2,3)
% plot(t,x(:,3));
% xlabel('Time [s]')
% ylabel('Angular velocity [rad/s]')
% title('Joint 1 angular velocity')
% hold on
% plot(t,2*cos(2*t),'k:','linewidth',2)
% legend('PD controler output','Desired anglular velocity')
% 
% subplot(2,2,4)
% plot(t,x(:,4));
% xlabel('Time [s]')
% ylabel('Angular velocity [rad/s]')
% title('Joint 2 angular velocity')
% hold on
% plot(t,cos(t),'k:','linewidth',2)
% legend('PD controler output','Desired anglular velocity')
%%
figure(1)
subplot(2,2,3)
T1=I1*x(:,3)-I2*x(:,4);
T2=I2*x(:,4);
plot(t,T1)
ylabel('Torque[N.m]')
xlabel('Time[s]')
title('Joint 1 adaptive control method')
subplot(2,2,4)
plot(t,T2)
ylabel('Torque[N.m]')
xlabel('Time[s]')
title('Joint 2 adaptive control method')