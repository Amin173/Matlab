clc
clear
U=eye(3);
Dt=[1;1;1];
mew0=[0;0;0];
sigma=1*eye(1);
mew=mew0;
for i=1:size(U,1)
    ut=U(i,:)';
    dt=Dt(i);
    [z x]=measurement_2(mew0,U,Dt,i);
    [mew,sigma] = EKF_underwater_robot(  mew,sigma,ut,z,dt);
figure()
hold on
error_ellipse(sigma,mew,'conf',0.68,'style','');
xlabel('x[m]')
ylabel('y[m]')
zlabel('z[m]')
title(['time =' num2str(sum(Dt(1:i))) '[s]']);
end
