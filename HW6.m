clear
clc
% y0=[-45*pi/180;45*pi/180;0;0;0;0;0;0;0];
y0=[-45*pi/180;45*pi/180;0;0;0.001;0.001;0.001;0.001;0.001];
tf=10;
global kd lambda T TT I1 I2

lambda=10;
kd=[300 0;0 100];
[t,q]=ode45(@equ2,[0 tf],y0);


% A=xy(q(:,1),q(:,2));
% hold on
% plot(A(:,1),A(:,2));
% xlabel('X [m]')
% ylabel('Y [m]')
% title('End effector trajectory')
% txt1 = '\leftarrow Actual start';
% txt2='\leftarrow Actual end point';
% text(A(1,1),A(1,2),txt1)
% text(A(end,1),A(end,2),txt2)
% a=2*sin(2*t);
% b=sin(t);
% B=xy(a,b);
% hold on
% plot(B(:,1),B(:,2),'k:','linewidth',2);
% legend('Followed Path','Desired path')
% txt1 = '\leftarrow Desired start';
% txt2='\leftarrow Desired end point';
% text(B(1,1),B(1,2),txt1)
% text(B(end,1),B(end,2),txt2)
% set(gca,'Fontsize',14)


% subplot(2,2,1)
% plot(t,q(:,1));
% xlabel('Time [s]')
% ylabel('Joint Angle [rad]')
% title('Joint 1 angular position')
% hold on
% plot(t,a,'k:','linewidth',2)
% legend('Adaptive controler output','Desired angle')
% 
% subplot(2,2,2)
% plot(t,q(:,2));
% xlabel('Time [s]')
% ylabel('Joint Angle [rad]')
% title('Joint 2 angular position')
% hold on
% plot(t,b,'k:','linewidth',2)
% legend('Adaptive controler output','Desired angle')
% 
% subplot(2,2,3)
% plot(t,q(:,3));
% xlabel('Time [s]')
% ylabel('Angular velocity [rad/s]')
% title('Joint 1 angular velocity')
% hold on
% plot(t,2*cos(2*t),'k:','linewidth',2)
% legend('Adaptive controler output','Desired anglular velocity')
% 
% subplot(2,2,4)
% plot(t,q(:,4));
% xlabel('Time [s]')
% ylabel('Angular velocity [rad/s]')
% title('Joint 2 angular velocity')
% hold on
% plot(t,cos(t),'k:','linewidth',2)
% legend('Adaptive controler output','Desired anglular velocity')
%% Joint Torques
figure(1)
subplot(2,2,1)
T1=I1*q(:,3)-I2*q(:,4);
T2=I2*q(:,4);
plot(t,T1)
ylabel('Torque[N.m]')
xlabel('Time[s]')
title('Joint 1 adaptive control method')
subplot(2,2,2)
plot(t,T2)
ylabel('Torque[N.m]')
xlabel('Time[s]')
title('Joint 2 adaptive control method')