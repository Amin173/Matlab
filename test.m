clc
clear
%%
% xData = [2,3,4,6];
% yData = [0,23,3,40]
A = [1, 2, 2^2, 2^3;
    1, 3, 3^2, 3^3;
    1, 4, 4^2, 4^3;
    1, 6, 6^2, 6^3];
b = [0; 23; 3; 40];
coef = A\b;
coef = flipud(coef);
x1 = [2,3,4,6];
y = [0,23,3,40]
xData = 2:0.1:6;
yData = polyval(coef,xData);

xVector = xData'; yVector = yData'; vectorOfOnes = ones(length(xData),1); 

A = [xVector vectorOfOnes]; b = yVector;
slopeAndOffset = (A'*A)\(A'*b); 
m = slopeAndOffset(1); y0 = slopeAndOffset(2);
scatter(xData,yData), hold on, x = 1.5:0.1:3; 
plot(x,m*x+y0) 