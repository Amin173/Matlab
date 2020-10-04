function [dydt] = equ3(t,q)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
k=4;
b=0.8;
vr=1;
wr=10*pi/180;
xr=[10;2;45*pi/180];
% xr=[-10;10;135*pi/180];
e=[cos(q(3)) sin(q(3)) 0;-sin(q(3)) cos(q(3)) 0;0 0 1]*(xr-[q(1);q(2);q(3)]);
a=(wr^2+b*vr^2)^0.5;
k1=2*k*a;
k4=b*abs(vr);
k3=2*k*a;
u=[-k1*e(1);-k4*vr*sin(e(3))*e(2)-k3*e(3)];
% u=[-k1*e(1);-k4*sign(vr)*e(2)-k3*e(3)];
dydt=[(vr*cos(e(3))-u(1))*cos(q(3));(vr*cos(e(3))-u(1))*sin(q(3));wr-u(2)];
end

