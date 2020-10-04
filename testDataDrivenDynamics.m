
nstate = 20;
ninputs = 1;
noutputs = 1

seed = 353423
rng(seed)
% generate random linear dynamical system
sys = drss(nstate,ninputs,noutputs);

A = sys.a;

%% impulse response
[y,t,x] = impulse(sys);

%%
X = x'; % transpose data to match convention

% identify model using least-squares fit to data
% this is the same as X(:,2:end)*pinv(X(:,1:end-1))
Adata = X(:,2:end)/X(:,1:end-1);

figure
subplot(1,2,1)
contourf(A)
title('True A')
subplot(1,2,2)
contourf(Adata)
title('A from data')

%%
figure
plot(eig(A),'kx')
hold on
plot(eig(Adata),'ro')
legend('True Eigenvalues','Predicted')
xlabel('Real(\lambda)')
ylabel('Imag(\lambda)')

%% predict time evolution of training data
Xpred = zeros(size(X));
Xpred(:,1:2) = X(:,1:2);
for ii = 2:(size(X,2)-1)
    Xpred(:,ii+1) = Adata*Xpred(:,ii);
end

ypred = sys.C*Xpred;

figure
plot(t,y)
hold on
plot(t,ypred,'--')
legend('True','Predicted')
xlabel('Time')
ylabel('y')

%%
figure
subplot(1,2,1)
contourf(X)
xlabel('Time')
ylabel('x')
title('True')
subplot(1,2,2)
contourf(Xpred)
xlabel('Time')
ylabel('x')
title('Predicted')

%% Predict response to a different initial condition

x0 = randn(size(sys.b));
[y,t,xtrue] = initial(sys,x0);
Xtrue = xtrue';

Xpred = zeros(size(Xtrue));
Xpred(:,1) = Xtrue(:,1);

for ii = 1:(size(Xtrue,2)-1)
    Xpred(:,ii+1) = Adata*Xpred(:,ii);
end

ypred = sys.C*Xpred;

figure
plot(t,y)
hold on
plot(t,ypred,'--')
legend('True','Predicted')
xlabel('Time')
ylabel('y')
