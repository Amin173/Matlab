clear
clc

y0=[0;0;0];
tf=10;
[t,q]=ode45(@equ4,[0 tf],y0);
x=q(:,1);
y=q(:,2);
teta=q(:,3);
figure(1)

hold on
subplot(1,3,1)
plot(t,x-10,'-o')
xlabel('t [s]')
ylabel('Time response of pose error [m]')
title('X error')
hold on
subplot(1,3,2)
plot(t,y-2,'-o')
title('Y error')
xlabel('t [s]')
ylabel('Time response of pose error [m]')
hold on
subplot(1,3,3)
hold on
plot(t,(180/pi)*teta-45,'-o')
title('Teta error')
xlabel('t [s]')
ylabel('Time response of pose error [deg]')
% figure()
% quiver(x,y,cos(teta),sin(teta))
% xlabel('X [m]')
% ylabel('Y [m]')
% axis equal
%%
clear
clc

y0=[0;1;10*pi/180];
tf=20;
[t,q]=ode45(@equ5,[0 tf],y0);
x=q(:,1);
y=q(:,2);
% figure()
% plot(t,x,t,y)
figure()
quiver(x,y,cos(q(:,3)),sin(q(:,3)));
axis equal
xlabel('X [m]')
ylabel('Y [m]')