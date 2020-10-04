clc 
clear
% r=10e-6:10e-6:100e-6;
r=70e-6;
Ap=(0.0254^2)*(1.7*2.5);
row=2500;
% k=0:0.1:1;
k=0.01;
A=k*Ap;
for i=1:length(r)
n=A./(4*pi*r(i).^2);
m=1000*n*row*(4/3)*pi*r(i).^3

% hold on
% plot(k,m);
end
% grid on
% grid minor
% xlabel('A/Ap');
% ylabel('m [gr]');
% legend('show')

