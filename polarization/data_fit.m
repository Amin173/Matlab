clc
clear
close all
%% Importing Data
time_dch_l=load('time_dch.mat');
time_dch=time_dch_l.time;

V_dch_l=load('V_exp_dch.mat');
V_dch=V_dch_l.V_exp;

V_ch_l=load('V_exp_ch.mat');
V_ch=V_ch_l.V_exp;

time_ch_l=load('time_ch.mat');
time_ch=time_ch_l.time;
%% Parameters Initialization
R=8.314;
T=298;
F=96487;
Q=R*T/F;
%% Data Interpolation
time=0:0.0001:5;
V_dch_inp=interp1(time_dch,V_dch,time,'linear');
V_ch_inp=interp1(time_ch,V_ch,time,'linear');

time_ch_inp=time(~isnan(V_ch_inp));
time_dch_inp=time(~isnan(V_dch_inp));
V_dch_inp=V_dch_inp(~isnan(V_dch_inp));
V_ch_inp=V_ch_inp(~isnan(V_ch_inp));
%% Least square non-linear data fitting
zp_dch=time_dch_inp;
zp_ch=time_ch_inp;
% The curve fitting parameters are
%x(1)=E0
%x(2)=A0
% x(3)=B0
% fun=@(x,zp)(x(1)+x(2)*exp(-x(3)*zp));
fun=@(x,zp)(x(1)*exp(x(2)*zp)+x(3)*exp(x(4)*zp));

lb=-inf*ones(1,4);
ub=inf*ones(1,4);

x0_dch=[1,1,1,1];
x0_ch=[1,1,1,1];

[x_dch,error_dch] = lsqcurvefit(fun,x0_dch,zp_dch,V_dch_inp);
[x_ch,error_ch] = lsqcurvefit(fun,x0_ch,zp_ch,V_ch_inp);

problem_ch = createOptimProblem('lsqcurvefit','x0',x0_ch,'objective',fun,...
        'lb',lb,'ub',ub,'xdata',zp_ch,'ydata',V_ch_inp);

    ms_ch = MultiStart('PlotFcns',@gsplotbestf);
[xmulti_ch,errormulti_ch] = run(ms_ch,problem_ch,50);

problem_dch = createOptimProblem('lsqcurvefit','x0',x0_dch,'objective',fun,...
        'lb',lb,'ub',ub,'xdata',zp_dch,'ydata',V_dch_inp);

    ms_dch = MultiStart('PlotFcns',@gsplotbestf);
[xmulti_dch,errormulti_dch] = run(ms_dch,problem_dch,50);

V_bnd_dch=fun(xmulti_dch,zp_dch);
V_bnd_ch=fun(xmulti_ch,zp_ch);

%% Plotting results
% times = linspace(zp(1),zp(end));
figure()
hold on
plot(time_ch,V_ch,'g:','LineWidth',2)
plot(zp_ch,V_bnd_ch,'k','LineWidth',2)

plot(time_dch,V_dch,'b:','LineWidth',2)
plot(zp_dch,V_bnd_dch,'r','LineWidth',2)



% for i=1:length(zd)
% line([zrev(i) zd(i)],[V_ch_rev(i) V_ch_rev(i)])
% % end
% 
% plot(z_rev',V_scaning','LineWidth',2)
xlabel('time [s]')
ylabel('V_exp [V]')
legend('Experimental charging','Fitted charging bnd curve','Experimental discharge','Fitted discharge bnd curve')
% title('Data and Fitted Curve')
%% save data
save('x','xmulti_ch','xmulti_dch')