
function [mew,sigma] = EKF_underwater_robot(  mew0,sigma0,ut,z,dt)

a1=0.001;
a2=0.001;
a3=0.001; 
x=mew0(1);
y=mew0(2);
z=mew0(3);
xd=ut(1);
yd=ut(2);
zd=ut(3);

Gt=eye(3);
Vt=eye(3)*dt;

Mt=[a1*xd^2,0,0;
    0,a2*yd^2,0;
    0,0,a3*zd^2];

mew_bar=mew0+[xd;yd;zd]*dt;

sigma_bar=Gt*sigma0*Gt'+Vt*Mt*Vt';
    

Q=(0.01)^2*eye(3);

mx=[10 10 15];
my=[10 -10 0];
mz=[10 20 -10];
% Beacons=[mx' my' mz'];

% for i=1:size(z,2)
for j=1:3
q(j)=norm([mx(j)-mew_bar(1) my(j)-mew_bar(2) mz(j)-mew_bar(3)]);
H(j,:)=[-(mx(j)-mew_bar(1))/q(j) -(my(j)-mew_bar(2))/q(j) -(mz(j)-mew_bar(3))/q(j)];
end
z_hat=q;

S=H*sigma_bar*H'+Q;
K=sigma_bar*H'*inv(S);
mew_bar=mew_bar+K*(z-z_hat)';
sigma_bar=(eye(3)-K*H)*sigma_bar;

% end

mew=mew_bar;
sigma=sigma_bar;

% x=[mew;sigma;pzt];
end

