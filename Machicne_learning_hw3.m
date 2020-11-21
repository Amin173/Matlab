
clc; clear; close all;

x = -20:0.1:20;
a1=0; b1=1; a2=1; b2=2;
f = exp(abs(x-a2)/b2 - abs(x-a1)/b1);
plot(x, f)
xlabel('x')
ylabel('Likelihood(x)')