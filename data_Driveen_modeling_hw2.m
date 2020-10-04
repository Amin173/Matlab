%% Problem 1.1
clc; clear;
y = zeros(40, 1);
b = load('hw2Q1_b.dat');
C = load('hw2Q1_C.dat');
% Solve y = Theta * s for "s"
n = size(C, 2); % dimension of s
p = length(b); % number of measurements, dim(y)
% L1 minimum norm solution s_L1
cvx_begin;
variable s_L1(n);
minimize( norm(s_L1,1) );
subject to
C*s_L1 == b;
cvx_end;
s_L1 = abs(s_L1)>1e-5;
ids = find(s_L1 == 1);
disp(['The fake coin indexes using l_1 minimzation are: {', num2str(ids'), '}'])
%% Problem 1.2
K = sum(s_L1);
tol = 1e-10;
maxiterations = 1e4;
s = cosaomp(C,b,K,tol,maxiterations);
s = abs(s)>1e-5;
ids = find(s == 1);
disp(['The fake coin indexes using cosaomp are: {', num2str(ids'), '}'])
%% Problem 2
clear all
close all
clc
dy  = 0.01;
y = (-2:dy:2)'; %spatial coordinate

dt = 0.1;
Nt = 101;
tend = dt*(Nt-1);
t = 0:dt:tend; %time

% define function
amp1 = 1;
y01 = 0.5;
sigmay1  = 0.6;

amp2 = 1.2;
y02 = -0.5;
sigmay2  = 0.3;
omega1 = 1.3;
omega2 = 4.1;

v1 = amp1*exp(-(y-y01).^2/(2*sigmay1^2));
v2 = amp2*exp(-(y-y02).^2/(2*sigmay2^2));
X = v1*exp(1i*omega1*t) + v2*exp(1i*omega2*t);
[U,S,V] = svd(X);
%% plot the right eigenvectors vs the spacial coordinates
hold on
plot(y, real(U(:,1)), y, real(U(:,2)))
xlabel('Spatial coordinates (y)')
ylabel('Right eigenvectos (U)')
legend('U1', 'U2')
%% compute the dynamic mode decomposition of X
r = sum(abs(S(:))>1e-5);
[Phi, Lambda, W, b] = DMD(X(:,1:end-1),X(:,2:end),r);
%%
plot(W(:,1))
plot(W(:,2))
%%
Lambda_c = log(Lambda)/dt;
%% compute the optimal sensor locations for two sensors
p = 2; % # of sensors p
Psi = U(:,1:r);
[Q,R,pivot] = qr(Psi','vector');
C = zeros(p,size(Psi,1));
for j=1:p
C(j,pivot(j))=1;
end

Theta = C*Psi;
ym = X(pivot(1:p),1); % Measure at pivot locations
a = Theta\ym; % Estimate coefficients
dataRecon = U(:,1:r)*a; % Reconstruct data
figure()
hold on
plot(real(dataRecon), 'o')
plot(X(:,1))
legend('Original data', 'Reconstructed data')
xlabel('Time[s]')
ylabel('Spacial coordinates')