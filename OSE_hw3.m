clc 
clear
gama=[0;1];
phi=[0.6 -0.2;0.2 0.6];
syms x1 x2 x12 W
simplify(phi*[x1 x12;x12 x2]*phi')
y=solve((9*x1)/25 + x2/25 - (6*x12)/25==x1,(3*x1)/25 - (3*x2)/25 + (8*x12)/25==x12, W+ x1/25 + (9*x2)/25 + (6*x12)/25==x2);
y.x1
y.x2
y.x12
[7/39 -3/13;-3/13 58/39]
X = dlyap(phi,gama*gama')