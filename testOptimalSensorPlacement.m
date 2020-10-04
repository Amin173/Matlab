% Comparison between compressed sensing and optimal sensor placement
% reconstruction on flow over a flat plate airfoil at
% an angle of attack of 30 degrees, at a Reynolds number of 100,
% using the method outlined in Manohar et al., "Data-Driven
% Sparse Sensor Placement for Reconstruction, IEEE Control Systems
% Magazine, 2018.

clear all
close all
clc

load DataVort30.mat
n = size(DataVortR,1);
Xrt = Xr';
Yrt = Yr';
Xvec = Xrt(:);
Yvec = Yrt(:);

%% Do SVD

% compute mean
DataVortMean = mean(DataVortR,2);
% subtract mean
DataVortSub = DataVortR - DataVortMean;

[U,S,V] = svd(DataVortSub,'econ');
%% Plot left singular vectors (POD modes)
figure
contourFactor = 0.5;
for ii = 1:8
    subplot(4,2,ii)
    makePlotNoLine(Xr,Yr,reshape(U(:,ii),size(Yr'))',contourFactor)
    %ylim([-2.5,2.5])
    plot(xBody,yBody,'k','Linewidth',2)
    title(['Mode ',num2str(ii)])
end

%% Try compressed sensing using random sensor locations
snapind = 1;
nSensors = 20;
nModesSparse = 20;
NoiseLevel = 0;
perm = round(rand(nSensors, 1) * n); % choose random sensor locations
y = DataVortSub(perm,snapind); % compressed measurement
y = y + NoiseLevel*randn(size(y)); % add noise, if desired
Theta = U(perm,1:nModesSparse);
s = cosamp(Theta,y,10,1.e-10,10);

%%
figure
contourFactor = 0.1;
subplot(2,1,1)
makePlotNoLine(Xr,Yr,reshape(DataVortMean+DataVortSub(:,snapind),size(Yr'))',contourFactor)
plot(xBody,yBody,'k','Linewidth',2)
plot(Xvec(perm),Yvec(perm),'rx')
title('True')
subplot(2,1,2)
makePlotNoLine(Xr,Yr,reshape(DataVortMean+U(:,1:nModesSparse)*s,size(Yr'))',contourFactor)
plot(xBody,yBody,'k','Linewidth',2)
%  plot(Xvec(perm),Yvec(perm),'rx')
title('Reconstructed, random sensors')

%% Do optimal sensor placement version
r = 10;

[~,~,pivot] = qr(U(:,1:r)','vector');
sensors = pivot(1:r);

% this is a = (C Phi_r)^-1 y:
y = DataVortSub(sensors,snapind);
y = y + NoiseLevel*randn(size(y));
% reconstruct state in low-dimensional space
a = U(sensors,1:r)\y;

% Estimate of full state
%this is x = Phi_r a
DataReconLS = U(:,1:r)*a;

%%
figure
contourFactor = 0.1;
subplot(2,1,1)
makePlotNoLine(Xr,Yr,reshape(DataVortMean+DataVortSub(:,snapind),size(Yr'))',contourFactor)
plot(xBody,yBody,'k','Linewidth',2)
plot(Xvec(sensors),Yvec(sensors),'kx')
title('True')
subplot(2,1,2)
makePlotNoLine(Xr,Yr,reshape(DataVortMean+DataReconLS,size(Yr'))',contourFactor)
plot(xBody,yBody,'k','Linewidth',2)
%  plot(Xvec(perm),Yvec(perm),'rx')
title('Reconstructed, optimal sensors')
%% Compare condition numbers

cond(Theta) % for random placement of sensors
cond(U(sensors,1:r)) % for optimal placement




%% try again (works also);
%{
C = QRsensors(U(:,1:r),r);
DataReconLS2 = U(:,1:r)*((C*U(:,1:r))\DataVortSub(sensors,snapind));

figure
contourFactor = 0.1;
    subplot(2,1,1)
    makePlotNoLine(Xr,Yr,reshape(DataVortMean+DataVortSub(:,snapind),size(Yr'))',contourFactor)
    plot(xBody,yBody,'k','Linewidth',2)
    plot(C*Xvec,C*Yvec,'rx')
    title('True')
    subplot(2,1,2)
    makePlotNoLine(Xr,Yr,reshape(DataVortMean+DataReconLS2,size(Yr'))',contourFactor)
    plot(xBody,yBody,'k','Linewidth',2)
  %  plot(Xvec(perm),Yvec(perm),'rx')
    title('Reconstructed, random sensors')
%}