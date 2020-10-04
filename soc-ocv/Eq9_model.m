clc
clear
close all
%% Importing Data
SOC_dch_l=load('SOC_dch.mat');
SOC_dch=SOC_dch_l.SOC_dch;

V_dch_l=load('V_dch.mat');
V_dch=V_dch_l.V_dch;

V_ch_l=load('V_ch.mat');
V_ch=V_ch_l.V_ch;

SOC_ch_l=load('SOC_ch.mat');
SOC_ch=SOC_ch_l.SOC_ch;

scaning02_l=load('scaning_curve_zrev_02');
SOC_scaning02=scaning02_l.SOC_ch;
V_scaning02=scaning02_l.V_ch;

%% Parameters Initialization
R=8.314;
T=298;
F=96487;
Q=R*T/F;
%% Data Interpolation
SOC=0.001:0.001:0.988;
V_dch_inp=interp1(SOC_dch,V_dch,SOC,'linear');
V_ch_inp=interp1(SOC_ch,V_ch,SOC,'linear');

SOC_ch_inp=SOC(~isnan(V_ch_inp));
SOC_dch_inp=SOC(~isnan(V_dch_inp));
V_dch_inp=V_dch_inp(~isnan(V_dch_inp));
V_ch_inp=V_ch_inp(~isnan(V_ch_inp));
%% Least square non-linear data fitting
zp_dch=SOC_dch_inp;
zp_ch=SOC_ch_inp;

% The curve fitting parameters are:
%x(1):E0
% x(2):n
%x(3):A0
% x(4):B0

%% Fitness function
% fun=@(x,zp)(x(1)+(Q)*log((zp)./(1-zp))+(0.5*Q)*(2*x(2)*(1-2*zp)+x(3)*(1-3*zp.^2)));
fun=@(x,zp)(x(1)+(Q/x(2))*log((zp)./(1-zp+1e-8))+0.5*Q*(-2*x(3)+4*x(3)*(1-zp)-2*x(4)+6*x(4)*(1-zp)-3*x(4)*(1-zp).^2));

%%

x0_dch=[3.278,0.8621,-0.1576,-0.04263];
x0_ch=[3.022,0.2168,8.057,27.33];


lb=-inf*ones(1,4);
ub=inf*ones(1,4);

[x_dch,error_dch] = lsqcurvefit(fun,x0_dch,zp_dch,V_dch_inp);
[x_ch,error_ch] = lsqcurvefit(fun,x0_ch,zp_ch,V_ch_inp);




%% zrev and zd from experimental data
% zrev=[0.8,0.6,0.4,0.2];
% 
% V_ch_rev=interp1(SOC_ch,V_ch,zrev,'linear');
% 
% [V, index] = unique(V_dch); 
% 
% zd=interp1(V,SOC_dch(index),V_ch_rev,'linear');
%% zrev and zd from fitted data
zrev=[0.8,0.6,0.4,0.2];

V_ch_rev=fun(x_ch,zrev);

[V, index] = unique(fun(x_dch,zp_dch)); 

zd=interp1(V,zp_dch(index),V_ch_rev,'linear');

%% Scanning curves coefficients

for i=1:length(zd)
    z_scaning(i,:)=linspace(0,zd(i),1000);
    z_rev(i,:)=linspace(0,zrev(i),1000);
end
% V_scaning=interp1(SOC_dch,V_dch,z_scaning,'linear');
V_scaning=fun(x_dch,z_scaning);
zp=0:0.001:0.999;
V_bnd_dch=fun(x_dch,zp);
V_bnd_ch=fun(x_ch,zp);

%% Plotting results
% times = linspace(zp(1),zp(end));
figure(1)
hold on
plot(SOC_dch,V_dch,'c:','LineWidth',2)
plot(SOC_ch,V_ch,'g:','LineWidth',2)
plot(zp,V_bnd_dch,'c','LineWidth',2)
plot(zp,V_bnd_ch,'g','LineWidth',2)

% for i=1:length(zd)
% line([zrev(i) zd(i)],[V_ch_rev(i) V_ch_rev(i)])
% end

plot(z_rev',V_scaning','LineWidth',2)
plot(SOC_scaning02,V_scaning02,'r:','LineWidth',2)
xlabel('SOC')
ylabel('OCV')
legend('Experimental discharge','Experimental charging','Fitted discharge bnd curve','Fitted charging bnd curve','Scaning curve No. 1','Scaning curve No. 2','Scaning curve No. 3','Scaning curve No. 4','Scaning Curve No. 4 Exp. data')
% title('Data and Fitted Curve')
%% Comparing experimental scaning curve with simulated scaning curve

V_dch_rev=fun(x_dch,zrev(4));

[V, index] = unique(fun(x_ch,zp_ch)); 

zd_ch=interp1(V,zp_ch(index),V_dch_rev,'linear');

%% Scanning curves coefficients

for i=1:length(zd)
    z_scaning_ch=linspace(0,zd_ch,1000);
    z_rev_ch=linspace(0,zrev(4),1000);
end
% V_scaning=interp1(SOC_dch,V_dch,z_scaning,'linear');
V_scaning_ch=fun(x_ch,z_rev_ch);

figure(2)
hold on
plot(z_rev(4,:)',V_scaning(4,:)','LineWidth',2)
plot(z_rev_ch,V_scaning_ch,'c','LineWidth',2);
plot(SOC_scaning02,V_scaning02,'r:','LineWidth',2)
xlabel('SOC')
ylabel('OCV')
legend('Scaning curve No. 4','Fitted charging bnd curve','Scaning Curve No. 4 Expe. data')