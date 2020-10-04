clc
clear
mew0=[80;100;0];
sigma0=[25^2 0 0;0 50^2 0;0 0 0];
ut=[10;5*pi/180];
dt=1;
t=9;
mew=mew0;
sigma=sigma0;
for i=1:dt:t

[z x]=measurement_1(t);
[mew sigma VMV GSG Q HSH S z_hat  mew_bar sigma_bar]=EKF_locolization_known_correspondences(mew,sigma,ut,z,dt);

end
z_hat(2)=z_hat(2)*180/pi;
%% Prediction Step of EKF
hold on
error_ellipse(sigma0(1:2,1:2),mew0(1:2),'conf',0.34,'style','');
error_ellipse(VMV(1:2,1:2),mew_bar(1:2),'conf',0.34,'style','');
error_ellipse(GSG(1:2,1:2),mew_bar(1:2),'conf',0.34,'style','');
error_ellipse(sigma_bar(1:2,1:2),mew_bar(1:2),'conf',0.34,'style','');

legend('Sigma(t-1)','V*M*V^T','G*Sigma(t-1)*G^T','Sigma(t)')
title('Prediction step of EKF algorithm - a1=0.1, a4=0.1')
xlabel('x[cm]')
ylabel('y[cm]')
axis equal
%% Measurement Prediction

hold on
title('Measurement Prediction - a1=0.1, a4=0.1')
subplot(1,2,1)
error_ellipse(sigma_bar(1:2,1:2),mew_bar(1:2),'conf',0.34,'style','');
xlabel('x[cm]')
ylabel('y[cm]')
legend('SigmaBar')
subplot(1,2,2)
hold on
error_ellipse(Q(1:2,1:2),z_hat(1:2),'conf',0.34,'style','');subplot(1,2,2)
error_ellipse(HSH(1:2,1:2),z_hat(1:2),'conf',0.34,'style','');
error_ellipse(S(1:2,1:2),z_hat(1:2),'conf',0.34,'style','');
legend('Q_t','H*SigmaBar*H^T','S_t')
xlabel('r[cm]')
ylabel('phi[deg]')
axis equal
%% Correction Step
hold on
error_ellipse(sigma_bar(1:2,1:2),mew_bar(1:2),'conf',0.34,'style','');
error_ellipse(sigma(1:2,1:2),mew(1:2),'conf',0.34,'style','');

legend('Sigma_bar(t-1)','Sigma')
title('Correction of EKF algorithm - a1=0.1, a4=0.1')
axis equal
xlabel('x[cm]')
ylabel('y[cm]')