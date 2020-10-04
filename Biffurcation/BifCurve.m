function f=BifCurve(cp,p,B)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
K=cp-B;
U=p-B;
N=100;
t=linspace(0,1,N);
r=norm(U)*t+norm(K)*(1-t);
T=cross(U,K);
eT=T/norm(T);
eU=U/norm(U);
% eK=K/norm(K);
alpha=atan2(norm(cross(U,T)),dot(U,T));
alpha=alpha*(alpha<=pi)+(2*pi-alpha)*(alpha>pi);
theta=4*alpha/N;
eTU=cross(eT,eU);
eTU=eTU/norm(eTU);
j=1:floor(N/4);
f=(r(4*j).*cos(j*theta)).*eU'+(r(4*j).*sin(j*theta)).*eTU'+B';
end

