function [dydt]=inversion(t,q)
M = 10;
m = 3;
b = 0.5;
l = 0.3;
I = (m*(2*l)^2)/12;
g = 9.8;
kp=100;
kd=30;
thd=pi;
wd=0;
p = I*(M+m)+M*m*l^2;
A = [-(I+m*l^2)*b/p  (m^2*g*l^2)/p ;
     -(m*l*b)/p       m*g*l*(M+m)/p ];
B = [(I+m*l^2)/p m*l/p;
     m*l/p (I+m*l^2)/p];
%     u=[-kp*(q(1)-xd)-kd(q(2)-xdotd);-kp*(q(3)-thd)-kd(q(4)-wd)];
    u=[-kp*(q(1)-0)-kd*(q(2)-0);-kp*(q(3)-thd)-kd*(q(4)-wd)];
V=A*[q(2);q(4)]+B*u;
dydt=[q(2);V(1);q(4);V(2)];


end