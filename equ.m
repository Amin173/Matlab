function dydt = equ(t,x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global kd kp lc2 l1 l2 lc1 I1 I2
m1=5;
m2=5;
l1=1;
l2=1;
g=0;
lc1=l1/2;
lc2=l2/2;
% T=[-2*exp(-3*t);-10*exp(-t)];
% T=[0;0];
T=[-kp*(x(1)-2*sin(2*t))-kd*(x(3)-4*cos(2*t));-kp*(x(2)-sin(t))-kd*(x(4)-cos(t))];
% J=[-l1*sin(x(1))-lc2*sin(x(1)+x(2)) -lc2*sin(x(1)+x(2));l1*cos(x(1))+lc2*cos(x(1)+x(2)) lc2*cos(x(1)+x(2))];
% T=-J.'*(kp*(xy(x(1),x(2))-xy(45*pi/180,10*pi/180))'+kd*J*[x(3);x(4)]);
I1=(1/12)*m1*l1^2;
I2=(1/12)*m2*l2^2;
H11=m1*lc1^2+I1+m2*(l1^2+lc2^2+2*l1*lc2*cos(x(2)))+I2;
H22=m2*lc2^2+I2;
H12=m2*l1*lc2*cos(x(2))+m2*lc2^2+I2;
h=m2*l1*lc2*sin(x(2));
G1=m1*lc1*g*cos(x(1))+m2*lc2*g*cos(x(1)+x(2))+m2*l1*g*cos(x(1));
G2=m2*lc2*g*cos(x(1)+x(1));
C=[-h*x(4) -h*(x(3)+x(4));h*x(3) 0];
H=[H11 H12;H12 H22];
V=H\(T-C*[x(3);x(4)]-[G1;G2]);
dydt=[x(3);x(4);V(1);V(2)];

end

