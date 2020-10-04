clc 
clear 
close all
h=3*2.25e-2; %height of the robot
r=2*2.25e-2; % radius of the robot
V=pi*r^2*h; % volume of a single robot
n=10; %number of robots per one vacuum pump
pf=0.7; % packing factor (particle volume ratio)
V_vacuum=n*pf*V; % air volume to be vacuumed [l]
flow_rate=0.1; %pump's flow rate [l/sec]
t=V_vacuum/flow_rate % estimated vacuum time

