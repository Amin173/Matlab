function [z x] = measurement_1( t )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
x0=80;
y0=100;
th0=0;
x=[x0;y0;th0];
dt=t/1000;
vt=10;
wt=5*pi/180;
for i=0:dt:t
theta=x(3);
x=x+[(vt/wt)*(-sin(theta)+sin(theta+wt*dt));(vt/wt)*(cos(theta)-cos(theta+wt*dt));wt*dt];
end
r=norm([x(1) x(2)]);
phi=atan2(x(2),x(1));
s=1;
z=[r;phi;s];
end

