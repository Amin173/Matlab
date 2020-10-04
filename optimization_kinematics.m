clc
clear
n = 12;
% x0 = 0:(2*pi/n);2*pi;
% dt_k = [0 0 0 0 0 0 pi pi pi pi pi pi pi];
% A = [];
% b = [];
% lb = zeros(n, 1);
% ub = ones(n, 1);
% Aeq = [cos(dt_k);sin(dt_k)];
% beq = 0;
% x = sym('x', [1 n]);
% dt_k = sym('dt_k', [1 n]);
syms l
x = [sym('x[0]'), sym('x[1]'), sym('x[2]'), sym('x[3]'), sym('x[4]'), sym('x[5]'), sym('x[6]'), sym('x[7]'), ...
    sym('x[8]'), sym('x[9]'), sym('x[10]'), sym('x[11]')];
dt_k = [sym('dt_k[0]'), sym('dt_k[1]'), sym('dt_k[2]'), sym('dt_k[3]'), sym('dt_k[4]'), sym('dt_k[5]'), ...
     sym('dt_k[6]'), sym('dt_k[7]'), sym('dt_k[8]'), sym('dt_k[9]'), sym('dt_k[10]'), sym('dt_k[11]')];

fun = (sum(((x(2:(end-1))-(x(3:end)+x(1:(end-2)))/2)/l-dt_k(2:end-1)).^2) +  (x(end)-((x(end-1)+x(1))/(2*l))-dt_k(end))^2 + (x(1)-((x(2)+x(end))/(2*l))-dt_k(1))^2);
for i= 1:n
jac(i,1) = diff(fun,x(i));
for j= 1:n
hess(i,j) = diff(jac(i), x(j));
end
end
% res = fmincon(fun,x0,A,b,Aeq,beq,lb,ub)
%% print jac
clc
for i = 1:n
disp([char(jac(i)), ','])
end
%% print hess
clc
for i =1:n
    fprintf('[')
    for j =1:n-1
        fprintf([char(hess(i,j)),','])
    end
    fprintf([char(hess(i,n)),'],\n'])
end