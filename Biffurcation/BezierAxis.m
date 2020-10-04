clc
clear
close all
%% Arbtrary input data
n=100; %n: number of data points in each segment (for simplicity, n is
% considered identical for each curve)

[x,y,z]=genData(5,n);
x1 = x(:,1);
x2 = x(:,2);
x3 = x(:,5);
x4 = x(:,4);
x5 = x(:,3);

y1 = y(:,1);
y2 = y(:,2);
y3 = y(:,5);
y4 = y(:,4);
y5 = y(:,3);

z1 = z(:,1);
z2 = z(:,2);
z3 = z(:,5);
z4 = z(:,4);
z5 = z(:,3);
