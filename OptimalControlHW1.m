clc; clear;
format shortEng
format compact
%% Grid definition
dx = 1;
dv = 1;
x = 0:dx:10;
v = 0:dv:5;
%% Parameter initializatiton
CD=0.4;
A= 0.1;
m=2;
row=1.204;
mew=0.2;
g=10;
q = 2; %weighting parameter for travel time cost funtion
M=1e9; %big number for cost penalty
%% Boundary conditions
J = zeros(length(v), length(x));
J(2:end, end) = M;
%% Initial conditions
Init = zeros(length(v), 1);
Init(2:end) = M; 
%% Cost function 
[X,Y] = meshgrid(v,v);
F = m*(abs(X.^2-Y.^2)/(2*dx))+0.5*CD*row*A*X.^2+mew*m*g; %energy cost function
JS = q * (2*dx)./(X.^2+Y.^2) + F; %total stage cost
a = M*((abs(X.^2-Y.^2)/(2*dx))>3); %admissible accelerations
%% Solution
for i= 0:(length(x)-3)
[JX, JY] = meshgrid(zeros(length(v),1), J(:, end-i));
J(:, end-i-1) = min(JS+JX+JY+a);
end
[JX, JY] = meshgrid(Init, J(:, 2));
J(:, 1) = min(JS+JX+JY+a);
%% Visualization
disp(J)
A = J==min(J);
vc=zeros(length(x),1);
for i=1:length(x)
    vc(i) = v(J(:,i)==min(J(:,i)));
end
ac=zeros(length(x),1);
for i=1:length(vc)-1
    ac(i)= (vc(i+1)^2-vc(i)^2)/dx;
end
ac(end)= (0-vc(end)^2)/dx;
hold on
plot(x,vc, '-s')
xlim([0 10])
% ylim([-3.5 3.5])