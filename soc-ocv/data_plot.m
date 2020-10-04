clc
clear
close all
%% Initializing variables
V_ch_l=load('V_ch.mat');
V_ch=V_ch_l.V_ch;

V_dch_l=load('V_dch.mat');
V_dch=V_dch_l.V_dch;

SOC_ch_l=load('SOC_ch.mat');
SOC_ch=SOC_ch_l.SOC_ch;

SOC_dch_l=load('SOC_dch.mat');
SOC_dch=SOC_dch_l.SOC_dch;

%% Ploting SOC-OCV in one figure
hold on
plot(SOC_ch,V_ch,'b','LineWidth',2)
plot(SOC_dch,V_dch,'r','LineWidth',2)
xlabel('SOC')
legend('Vch','Vdch')
%% 

x=0:0.01:1;

V_ch_inp=interp1(SOC_ch,V_ch,x,'linear');
V_dch_inp=interp1(SOC_dch,V_dch,x,'linear');


%% Ploting SOC-OCV in one figure
hold on
plot(x,V_ch_inp,'b','LineWidth',2)
plot(x,V_dch_inp,'r','LineWidth',2)
xlabel('SOC')
ylabel('OCV')
legend('V_ch_fit','V_dch_fit')
%% Hystersis
offset=V_ch_inp-V_dch_inp;
figure()
plot(x,offset)
xlabel('SOC')
ylabel('Offset potential')