clc
clear
alpha=0.1325;
T=1/alpha^0.5;
num=5*[1 1/T];
den=conv([1 0 0],[1 1/(alpha*T)]);
bode(num,den);
grid