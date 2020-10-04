function [mew sigma VMV GSG Q HSH S z_hat m_b s_b] = EKF_locolization_known_correspondences(  mew0,sigma0,ut,z,dt)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
a1=0.1;
a2=0.05;
a3=0.05; 
a4=0.1;
x=mew0(1);
y=mew0(2);
theta=mew0(3);
vt=ut(1);
wt=ut(2);
Gt=[1 0 -(vt/wt)*(cos(theta)-cos(theta+wt*dt));
    0 1 -(vt/wt)*(sin(theta)-sin(theta+wt*dt));
    0 0 1];
Vt=[(-sin(theta)+sin(theta+wt*dt))/wt  (vt/wt)*((sin(theta)-sin(theta+wt*dt))/wt+ cos(theta+wt*dt)*dt);
    (cos(theta)-cos(theta+wt*dt))/wt   (vt/wt)*(-(cos(theta)-cos(theta+wt*dt))/wt+ sin(theta+wt*dt)*dt);
    0 dt];

Mt=[a1*vt^2+a2*wt^2,0;
    0,a3*vt^2+a4*wt^2];

mew_bar=mew0+[(vt/wt)*(-sin(theta)+sin(theta+wt*dt));(vt/wt)*(cos(theta)-cos(theta+wt*dt));wt*dt];
m_b=mew_bar;  %to save the un-corrected mew_bar 
VMV=Vt*Mt*Vt';
GSG=Gt*sigma0*Gt';

sigma_bar=GSG+VMV;  
s_b=sigma_bar;%to save the un-corrected sigma_bar
    

Q=(80^2)*eye(3);

mx=[300];
my=[450];
ms=[1];

for i=1:size(z,2)
q=norm([mx(i)-mew_bar(1) my(i)-mew_bar(2)]);
z_hat=[q;atan2(my(i)-mew_bar(2),mx(i)-mew_bar(1))-mew_bar(3);ms(i)];
H=[-(mx(i)-mew_bar(1))/q -(my(i)-mew_bar(2))/q 0;(my(i)-mew_bar(2))/q -(mx(i)-mew_bar(1))/q -1;0 0 0];
HSH=H*sigma_bar*H';
S=HSH+Q;
K=sigma_bar*H'*inv(S);
mew_bar=mew_bar+K*(z-z_hat);
sigma_bar=(eye(3)-K*H)*sigma_bar;
end

mew=mew_bar;
sigma=sigma_bar;

% x=[mew;sigma;pzt];
end

