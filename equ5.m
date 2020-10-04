function [dydt] = equ3(t,q)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
k=0.8;
a=1;
v=1;

cs=0;
sdot=v*cos(q(3))*(1/(1-cs*q(2)));
ldot=v*q(3);
k2=a^2;
k3=2*k*a;
tetatildedot=-k2*v*q(2)-k3*abs(v)*q(3);


dydt=[sdot;ldot;tetatildedot];
end