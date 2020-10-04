clear
clc

y0=[0;0;pi/4;0];
tf=20;
[t,q]=ode45(@inversion2,[0 tf],y0);
th=q(:,3);
w=q(:,4);




subplot(2,2,1)
plot(t,q(:,1));
xlabel('Time [s]')
ylabel('x [m]')
title('Position along x axis [m]')
hold on
plot(t,0*t./t,'k:','linewidth',2)
legend('PD controler output','Desired postion')

subplot(2,2,2)
plot(t,q(:,2));
xlabel('Time [s]')
ylabel('V [m/s]')
title('Velocity in x direction [m/s]')
hold on
plot(t,0*t./t,'k:','linewidth',2)
legend('PD controler output','Desired velocity')

subplot(2,2,3)
plot(t,q(:,3));
xlabel('Time [s]')
ylabel('Theta [rad]')
title('Angular position Theta [rad]')
hold on
plot(t,pi*t./t,'k:','linewidth',2)
legend('PD controler output','Desired position')

subplot(2,2,4)
plot(t,q(:,4));
xlabel('Time [s]')
ylabel('Theta-dot [rad/s]')
title('Angular velocity Theta-dot [rad]')
hold on
plot(t,0*t./t,'k:','linewidth',2)
legend('PD controler output','Desired velocity')
