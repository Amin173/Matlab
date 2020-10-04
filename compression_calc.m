clc
clear
close all

l=0.1;
r=l/100:l/100:l/2;
foam_coeff=0.5;
comp_coeff=0.1:0.1:1;

P0=100;
Vcube=l^3;

for i=1:length(r)
n=(l/(2*r(i)))^2;
Vcyl(i)=n*pi*r(i)^2*l;
Vr=Vcube-Vcyl(i);
Va=foam_coeff*Vr;
Va2=Va+comp_coeff*Vcyl(i);
DP(:,i)=((Va./Va2)-1)*P0;

end
plot(comp_coeff,DP,'-*');
grid on
xlabel('Compression effieceny of the cylenders')
ylabel('Pressure differnce [kPa]')
% we don't have a limitation in meeting the requirments of volume variation

% %% Maximum motor load with 1cm shaft and 45deg thread
% % F=100e3*pi*0.025^2;
% % tow=F*0.5;
% % % result: Maximum pressure achievable is araound 20 [kPa]