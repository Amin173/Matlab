
clear all;
close all;
clc; 
t1=1.6*10^(-6);
T=20*10^(-6);
m=0.62;
Ap=12.5*(-(1-exp(-1333.333*(T-t1)))/(1-exp(-1333.33*t1-1333.33*(T-t1))));
Cp=12.5*(1+(exp(-1333.33*t1)*((exp(-1333.33*(T-t1))-1)/(1-exp(-1333.33*t1-1333.33*(T-t1))))));
An=12.5*((1-exp(-1333.33*(T-t1)))/(1-exp(-1333.33*t1-1333.33*(T-t1))));
Cn=12.5*((((1-exp(-1333.33*(T-t1)))*exp(-1333.33*t1))/(1-exp(-1333.33*t1-1333.33*(T-t1))))-1);
iav_p=-37.5*(Ap*exp(-1333.33*t1)-(Ap+Cp)+Cp*exp(-1333.33*(T-t1))-16666.625*t1);
iav_n=-37.5*(Ap*exp(-1333.33*t1)-(Ap+Cp)+Cp*exp(-1333.33*(T-t1))+16666.625*t1);
F_fk=1.2152;
% F_fk=1.8228;
% F_fk=2.4304;
% F_fk=3.038;
n=1;

for t=0:0.001:1
x(n)=9.25-9.25*cos(4*pi*t);
  if (x(n)>=0)&&(x(n)<0.9)
    R1(n)=sqrt((x(n)+5.2)^2+2.56);
    R2(n)=sqrt((x(n)-0.9)^2+2.56);
    R3(n)=sqrt((x(n)-7)^2+2.56);
    R4(n)=sqrt((x(n)-13.1)^2+2.56);
    alfa1(n)=acos(1.6/R1(n));
    alfa2(n)=acos(1.6/R2(n));
    alfa3(n)=acos(1.6/R3(n));
    alfa4(n)=acos(1.6/R4(n));
    g1(n)=inv(diag(diag((R1(n)')*(R1(n)))));
    g2(n)=inv(diag(diag((R2(n)')*(R2(n)))));
    g3(n)=inv(diag(diag((R3(n)')*(R3(n)))));
    g4(n)=inv(diag(diag((R4(n)')*(R4(n)))));
    B_total_1(n)=diag(2.4+((0.4*pi*iav_p)*g1(n)));
    B_total_2(n)=diag(2.4-((0.4*pi*iav_n)*g2(n)));
    B_total_3(n)=diag(2.4-((0.4*pi*iav_n)*g3(n)));
    B_total_4(n)=diag(2.4-((0.4*pi*iav_n)*g4(n)));
    k1(n)=diag(diag(B_total_1(n))/diag(B_total_2(n)));
    k3(n)=diag(diag(B_total_3(n))/diag(B_total_2(n)));
    k4(n)=diag(diag(B_total_4(n))/diag(B_total_2(n)));
    f1(n)=1*21.3376*iav_p*sin(alfa2(n))+k1(n)*21.3376*iav_p*sin(alfa1(n))+k3(n)*sin(alfa3(n))*21.3376*iav_p+k4(n)*sin(alfa4(n))*21.3376*iav_p-F_fk;
    delta_f1(n)=63.78-f1(n);
    ip1(n)=(delta_f1(n)+21.33*k1(n))/(21.33*k1(n));
    f1_new(n)=1*21.3376*iav_p*sin(alfa2(n))+k1(n)*21.3376*ip1(n)*sin(alfa1(n))+k3(n)*sin(alfa3(n))*21.3376*iav_p+k4(n)*sin(alfa4(n))*21.3376*iav_p-F_fk;
 
      elseif x(n)==0.9
    R1(n)=sqrt((x(n)+5.2)^2+2.56);
    R2(n)=sqrt((x(n)-0.9)^2+2.56);
    R3(n)=sqrt((x(n)-7)^2+2.56);
    R4(n)=sqrt((x(n)-13.1)^2+2.56);
    alfa1(n)=acos(1.6/R1(n));
    alfa2(n)=acos(1.6/R2(n));
    alfa3(n)=acos(1.6/R3(n));
    alfa4(n)=acos(1.6/R4(n));
    g1(n)=inv(diag(diag((R1(n)')*(R1(n)))));
    g2(n)=inv(diag(diag((R2(n)')*(R2(n)))));
    g3(n)=inv(diag(diag((R3(n)')*(R3(n)))));
    g4(n)=inv(diag(diag((R4(n)')*(R4(n)))));
    B_total_1(n)=diag(2.4+((0.4*pi*iav_p)*g1(n)));
    B_total_2(n)=2.4+eye(size((g2(n))));
    B_total_3(n)=diag(2.4-((0.4*pi*iav_n)*g3(n)));
    B_total_4(n)=diag(2.4-((0.4*pi*iav_n)*g4(n)));
    k1(n)=diag(diag(B_total_1(n))/diag(B_total_3(n)));
    k2(n)=diag(diag(B_total_2(n))/diag(B_total_3(n)));
    k4(n)=diag(diag(B_total_4(n))/diag(B_total_3(n)));
    f2(n)=1*21.3376*iav_p*sin(alfa3(n))+k1(n)*21.3376*iav_p*sin(alfa1(n))+k2(n)*sin(alfa2(n))*21.3376*iav_p+k4(n)*sin(alfa4(n))*21.3376*iav_p-F_fk;
  elseif (x(n)<=3.354)&&(x(n)>0.9)
    R1(n)=sqrt((x(n)+5.2)^2+2.56);
    R2(n)=sqrt((x(n)-0.9)^2+2.56);
    R3(n)=sqrt((x(n)-7)^2+2.56);
    R4(n)=sqrt((x(n)-13.1)^2+2.56);
    alfa1(n)=acos(1.6/R1(n));
    alfa2(n)=acos(1.6/R2(n));
    alfa3(n)=acos(1.6/R3(n));
    alfa4(n)=acos(1.6/R4(n));
    g1(n)=inv(diag(diag((R1(n)')*(R1(n)))));
    g2(n)=inv(diag(diag((R2(n)')*(R2(n)))));
    g3(n)=inv(diag(diag((R3(n)')*(R3(n)))));
    g4(n)=inv(diag(diag((R4(n)')*(R4(n)))));
    B_total_1(n)=diag(2.4+((0.4*pi*iav_p)*g1(n)));
    B_total_2(n)=diag(2.4+((0.4*pi*iav_p)*g2(n)));
    B_total_3(n)=diag(2.4-((0.4*pi*iav_n)*g3(n)));
    B_total_4(n)=diag(2.4-((0.4*pi*iav_n)*g4(n)));
    k2(n)=diag(diag(B_total_2(n))/diag(B_total_3(n)));
    k1(n)=diag(diag(B_total_1(n))/diag(B_total_3(n)));
    k4(n)=diag(diag(B_total_4(n))/diag(B_total_3(n)));
    f3(n)=1*21.3376*iav_p*sin(alfa3(n))+k1(n)*21.3376*iav_p*sin(alfa1(n))+k2(n)*sin(alfa2(n))*21.3376*iav_p+k4(n)*sin(alfa4(n))*21.3376*iav_p-F_fk;
    delta_f3(n)=77.62-f3(n);
    ip3_first(n)=(delta_f3(n)+21.33*k1(n))/(21.33*k1(n))-0.665;
    f3_new_first(n)=1*21.3376*iav_p*sin(alfa3(n))+k1(n)*21.3376*ip3_first(n)*sin(alfa1(n))+k2(n)*sin(alfa2(n))*21.3376*iav_p+k4(n)*sin(alfa4(n))*21.3376*iav_p-F_fk;
  elseif (x(n)<7)&&(x(n)>3.354)
    R1(n)=sqrt((x(n)+5.2)^2+2.56);
    R2(n)=sqrt((x(n)-0.9)^2+2.56);
    R3(n)=sqrt((x(n)-7)^2+2.56);
    R4(n)=sqrt((x(n)-13.1)^2+2.56);
    alfa1(n)=acos(1.6/R1(n));
    alfa2(n)=acos(1.6/R2(n));
    alfa3(n)=acos(1.6/R3(n));
    alfa4(n)=acos(1.6/R4(n));
    g1(n)=inv(diag(diag((R1(n)')*(R1(n)))));
    g2(n)=inv(diag(diag((R2(n)')*(R2(n)))));
    g3(n)=inv(diag(diag((R3(n)')*(R3(n)))));
    g4(n)=inv(diag(diag((R4(n)')*(R4(n)))));
    B_total_1(n)=diag(2.4+((0.4*pi*iav_p)*g1(n)));
    B_total_2(n)=diag(2.4+((0.4*pi*iav_p)*g2(n)));
    B_total_3(n)=diag(2.4-((0.4*pi*iav_n)*g3(n)));
    B_total_4(n)=diag(2.4-((0.4*pi*iav_n)*g4(n)));
    k2(n)=diag(diag(B_total_2(n))/diag(B_total_3(n)));
    k1(n)=diag(diag(B_total_1(n))/diag(B_total_3(n)));
    k4(n)=diag(diag(B_total_4(n))/diag(B_total_3(n)));
    f3(n)=1*21.3376*iav_p*sin(alfa3(n))+k1(n)*21.3376*iav_p*sin(alfa1(n))+k2(n)*sin(alfa2(n))*21.3376*iav_p+k4(n)*sin(alfa4(n))*21.3376*iav_p-F_fk;
    delta_f3(n)=77.62-f3(n);
    ip3_second(n)=(delta_f3(n)+21.33*k2(n))/(21.33*k2(n))-0.735;
    f3_new_second(n)=1*21.3376*iav_p*sin(alfa3(n))+k1(n)*21.3376*iav_p*sin(alfa1(n))+k2(n)*sin(alfa2(n))*21.3376*ip3_second(n)+k4(n)*sin(alfa4(n))*21.3376*iav_p-F_fk;
  elseif x(n)==7
    R1(n)=sqrt((x(n)+5.2)^2+2.56);
    R2(n)=sqrt((x(n)-0.9)^2+2.56);
    R3(n)=sqrt((x(n)-7)^2+2.56);
    R4(n)=sqrt((x(n)-13.1)^2+2.56);
    alfa1(n)=acos(1.6/R1(n));
    alfa2(n)=acos(1.6/R2(n));
    alfa3(n)=acos(1.6/R3(n));
    alfa4(n)=acos(1.6/R4(n));
    g1(n)=inv(diag(diag((R1(n)')*(R1(n)))));
    g2(n)=inv(diag(diag((R2(n)')*(R2(n)))));
    g3(n)=inv(diag(diag((R3(n)')*(R3(n)))));
    g4(n)=inv(diag(diag((R4(n)')*(R4(n)))));
    B_total_1(n)=diag(2.4+((0.4*pi*iav_p)*g1(n)));
    B_total_2(n)=diag(2.4+((0.4*pi*iav_p)*g2(n)));
    B_total_3(n)=2.4+eye(size((g3(n))));
    B_total_4(n)=diag(2.4-((0.4*pi*iav_n)*g4(n)));
    k2(n)=diag(diag(B_total_2(n))/diag(B_total_4(n)));
    k3(n)=diag(diag(B_total_3(n))/diag(B_total_4(n)));
    k1(n)=diag(diag(B_total_1(n))/diag(B_total_4(n)));
    f4(n)=1*21.3376*iav_p*sin(alfa4(n))+k1(n)*21.3376*iav_p*sin(alfa1(n))+k2(n)*sin(alfa2(n))*21.3376*iav_p+k3(n)*sin(alfa3(n))*21.3376*iav_p-F_fk;
  elseif (x(n)>7)&&(x(n)<9.25)
    R1(n)=sqrt((x(n)+5.2)^2+2.56);
    R2(n)=sqrt((x(n)-0.9)^2+2.56);
    R3(n)=sqrt((x(n)-7)^2+2.56);
    R4(n)=sqrt((x(n)-13.1)^2+2.56);
    alfa1(n)=acos(1.6/R1(n));
    alfa2(n)=acos(1.6/R2(n));
    alfa3(n)=acos(1.6/R3(n));
    alfa4(n)=acos(1.6/R4(n));
    g1(n)=inv(diag(diag((R1(n)')*(R1(n)))));
    g2(n)=inv(diag(diag((R2(n)')*(R2(n)))));
    g3(n)=inv(diag(diag((R3(n)')*(R3(n)))));
    g4(n)=inv(diag(diag((R4(n)')*(R4(n)))));
    B_total_1(n)=diag(2.4+((0.4*pi*iav_p)*g1(n)));
    B_total_2(n)=diag(2.4+((0.4*pi*iav_p)*g2(n)));
    B_total_3(n)=diag(2.4+((0.4*pi*iav_p)*g3(n)));
    B_total_4(n)=diag(2.4-((0.4*pi*iav_n)*g4(n)));
    k2(n)=diag(diag(B_total_2(n))/diag(B_total_4(n)));
    k3(n)=diag(diag(B_total_3(n))/diag(B_total_4(n)));
    k1(n)=diag(diag(B_total_1(n))/diag(B_total_4(n)));
    f5(n)=1*21.3376*iav_p*sin(alfa4(n))+k1(n)*21.3376*iav_p*sin(alfa1(n))+k2(n)*sin(alfa2(n))*21.3376*iav_p+k3(n)*sin(alfa3(n))*21.3376*iav_p-F_fk;
    delta_f5(n)=77.57-f5(n);
    ip5(n)=(delta_f5(n)+21.33*k2(n))/(21.33*k2(n))-0.630;
    f5_new(n)=1*21.3376*iav_p*sin(alfa4(n))+k1(n)*21.3376*iav_p*sin(alfa1(n))+k2(n)*sin(alfa2(n))*21.3376*ip5(n)+k3(n)*sin(alfa3(n))*21.3376*iav_p-F_fk;

   elseif x(n)==9.25
    f6(n)=0;
  elseif (x(n)<=12)&&(x(n)>9.25)
    R1(n)=sqrt((x(n)+5.2)^2+2.56);
    R2(n)=sqrt((x(n)-0.9)^2+2.56);
    R3(n)=sqrt((x(n)-7)^2+2.56);
    R4(n)=sqrt((x(n)-13.1)^2+2.56);
    alfa1(n)=acos(1.6/R1(n));
    alfa2(n)=acos(1.6/R2(n));
    alfa3(n)=acos(1.6/R3(n));
    alfa4(n)=acos(1.6/R4(n));
    g1(n)=inv(diag(diag((R1(n)')*(R1(n)))));
    g2(n)=inv(diag(diag((R2(n)')*(R2(n)))));
    g3(n)=inv(diag(diag((R3(n)')*(R3(n)))));
    g4(n)=inv(diag(diag((R4(n)')*(R4(n)))));
    B_total_1(n)=diag(2.4+((0.4*pi*iav_n)*g1(n)));
    B_total_2(n)=diag(2.4+((0.4*pi*iav_n)*g2(n)));
    B_total_3(n)=diag(2.4+((0.4*pi*iav_n)*g3(n)));
    B_total_4(n)=diag(2.4-((0.4*pi*iav_p)*g4(n)));
    k2(n)=diag(diag(B_total_2(n))/diag(B_total_4(n)));
    k3(n)=diag(diag(B_total_3(n))/diag(B_total_4(n)));
    k1(n)=diag(diag(B_total_1(n))/diag(B_total_4(n)));
    f7(n)=-(1*21.3376*iav_p*sin(alfa4(n))+k1(n)*21.3376*iav_p*sin(alfa1(n))+k2(n)*sin(alfa2(n))*21.3376*iav_p+k3(n)*sin(alfa3(n))*21.3376*iav_p-F_fk);
    delta_f7(n)=83.08+f7(n);
    ip7(n)=(delta_f7(n)+21.33*k3(n))/(21.33*k3(n));
    f7_new_first(n)=-(1*21.3376*(iav_p-0.3)*sin(alfa4(n))+k1(n)*21.3376*(iav_p-0.3)*sin(alfa1(n))+k2(n)*sin(alfa2(n))*21.3376*(iav_p-0.3)+k3(n)*sin(alfa3(n))*21.3376*ip7(n)-F_fk);
 elseif (x(n)<13.1)&&(x(n)>12)
    R1(n)=sqrt((x(n)+5.2)^2+2.56);
    R2(n)=sqrt((x(n)-0.9)^2+2.56);
    R3(n)=sqrt((x(n)-7)^2+2.56);
    R4(n)=sqrt((x(n)-13.1)^2+2.56);
    alfa1(n)=acos(1.6/R1(n));
    alfa2(n)=acos(1.6/R2(n));
    alfa3(n)=acos(1.6/R3(n));
    alfa4(n)=acos(1.6/R4(n));
    g1(n)=inv(diag(diag((R1(n)')*(R1(n)))));
    g2(n)=inv(diag(diag((R2(n)')*(R2(n)))));
    g3(n)=inv(diag(diag((R3(n)')*(R3(n)))));
    g4(n)=inv(diag(diag((R4(n)')*(R4(n)))));
    B_total_1(n)=diag(2.4+((0.4*pi*iav_n)*g1(n)));
    B_total_2(n)=diag(2.4+((0.4*pi*iav_n)*g2(n)));
    B_total_3(n)=diag(2.4+((0.4*pi*iav_n)*g3(n)));
    B_total_4(n)=diag(2.4-((0.4*pi*iav_p)*g4(n)));
    k2(n)=diag(diag(B_total_2(n))/diag(B_total_4(n)));
    k3(n)=diag(diag(B_total_3(n))/diag(B_total_4(n)));
    k1(n)=diag(diag(B_total_1(n))/diag(B_total_4(n)));
    f7(n)=-(1*21.3376*iav_p*sin(alfa4(n))+k1(n)*21.3376*iav_p*sin(alfa1(n))+k2(n)*sin(alfa2(n))*21.3376*iav_p+k3(n)*sin(alfa3(n))*21.3376*iav_p-F_fk);
    delta_f7(n)=83.08+f7(n);
    ip7(n)=(delta_f7(n)+21.33*k2(n))/(21.33*k2(n))-0.39;
    f7_new_second(n)=-(1*21.3376*(iav_p-0.15)*sin(alfa4(n))+k1(n)*21.3376*(iav_p-0.15)*sin(alfa1(n))+k2(n)*sin(alfa2(n))*21.3376*ip7(n)+k3(n)*sin(alfa3(n))*21.3376*(iav_p-0.15)-F_fk);

      elseif x(n)==13.1
    R1(n)=sqrt((x(n)+5.2)^2+2.56);
    R2(n)=sqrt((x(n)-0.9)^2+2.56);
    R3(n)=sqrt((x(n)-7)^2+2.56);
    R4(n)=sqrt((x(n)-13.1)^2+2.56);
    alfa1(n)=acos(1.6/R1(n));
    alfa2(n)=acos(1.6/R2(n));
    alfa3(n)=acos(1.6/R3(n));
    alfa4(n)=acos(1.6/R4(n));
    g1(n)=inv(diag(diag((R1(n)')*(R1(n)))));
    g2(n)=inv(diag(diag((R2(n)')*(R2(n)))));
    g3(n)=inv(diag(diag((R3(n)')*(R3(n)))));
    g4(n)=inv(diag(diag((R4(n)')*(R4(n)))));
    B_total_1(n)=diag(2.4+((0.4*pi*iav_n)*g1(n)));
    B_total_2(n)=diag(2.4+((0.4*pi*iav_n)*g2(n)));
    B_total_3(n)=diag(2.4+((0.4*pi*iav_n)*g3(n)));
    B_total_4(n)=2.4+eye(size((g4(n))));
    k2(n)=diag(diag(B_total_2(n))/diag(B_total_4(n)));
    k3(n)=diag(diag(B_total_3(n))/diag(B_total_4(n)));
    k1(n)=diag(diag(B_total_1(n))/diag(B_total_4(n)));
    f8(n)=-(1*21.3376*iav_p*sin(alfa4(n))+k1(n)*21.3376*iav_p*sin(alfa1(n))+k2(n)*sin(alfa2(n))*21.3376*iav_p+k3(n)*sin(alfa3(n))*21.3376*iav_p-F_fk);
  elseif (x(n)<=18.5)&&(x(n)>13.1)
    R1(n)=sqrt((x(n)+5.2)^2+2.56);
    R2(n)=sqrt((x(n)-0.9)^2+2.56);
    R3(n)=sqrt((x(n)-7)^2+2.56);
    R4(n)=sqrt((x(n)-13.1)^2+2.56);
    alfa1(n)=acos(1.6/R1(n));
    alfa2(n)=acos(1.6/R2(n));
    alfa3(n)=acos(1.6/R3(n));
    alfa4(n)=acos(1.6/R4(n));
    g1(n)=inv(diag(diag((R1(n)')*(R1(n)))));
    g2(n)=inv(diag(diag((R2(n)')*(R2(n)))));
    g3(n)=inv(diag(diag((R3(n)')*(R3(n)))));
    g4(n)=inv(diag(diag((R4(n)')*(R4(n)))));
    B_total_1(n)=diag(2.4+((0.4*pi*iav_n)*g1(n)));
    B_total_2(n)=diag(2.4+((0.4*pi*iav_n)*g2(n)));
    B_total_3(n)=diag(2.4+((0.4*pi*iav_n)*g3(n)));
    B_total_4(n)=diag(2.4+((0.4*pi*iav_n)*g4(n)));
    k2(n)=diag(diag(B_total_2(n))/diag(B_total_4(n)));
    k3(n)=diag(diag(B_total_3(n))/diag(B_total_4(n)));
    k1(n)=diag(diag(B_total_1(n))/diag(B_total_4(n)));
    f9(n)=-(1*21.3376*(iav_p)*sin(alfa4(n))+k1(n)*21.3376*(iav_p)*sin(alfa1(n))+k2(n)*sin(alfa2(n))*21.3376*(iav_p)+k3(n)*sin(alfa3(n))*21.3376*(iav_p)-F_fk);
    delta_f9(n)=65.37+f9(n);
    ip9(n)=(delta_f9(n)+21.33*k3(n))/(21.33*k3(n))-0.05;
    f9_new(n)=-(1*21.3376*(iav_p)*sin(alfa4(n))+k1(n)*21.3376*(iav_p)*sin(alfa1(n))+k2(n)*sin(alfa2(n))*21.3376*iav_p+k3(n)*sin(alfa3(n))*21.3376*(ip9(n))-F_fk);

  end
n=n+1;
end
f1_new(1,(length(f1_new))+1:(length(x)))=zeros(1,(length(x))-(length(f1_new)));
f3_new_first(1,(length(f3_new_first))+1:(length(x)))=zeros(1,(length(x))-(length(f3_new_first)));
f3_new_second(1,(length(f3_new_second))+1:(length(x)))=zeros(1,(length(x))-(length(f3_new_second)));
f5_new(1,(length(f5_new))+1:(length(x)))=zeros(1,(length(x))-(length(f5_new)));
f7_new_first(1,(length(f7_new_first))+1:(length(x)))=zeros(1,(length(x))-(length(f7_new_first)));
f7_new_second(1,(length(f7_new_second))+1:(length(x)))=zeros(1,(length(x))-(length(f7_new_second)));
f9_new(1,(length(f9_new))+1:(length(x)))=zeros(1,(length(x))-(length(f9_new)));


% f1(1,(length(f1))+1:(length(f9)))=zeros(1,(length(f9))-(length(f1)));
% f3(1,(length(f3))+1:(length(f9)))=zeros(1,(length(f9))-(length(f3)));
% f5(1,(length(f5))+1:(length(f9)))=zeros(1,(length(f9))-(length(f5)));
% f6(1,(length(f6))+1:(length(f9)))=zeros(1,(length(f9))-(length(f6)));
% f7(1,(length(f7))+1:(length(f9)))=zeros(1,(length(f9))-(length(f7)));
ft=f1_new+f3_new_first+f3_new_second+f5_new+f7_new_first+f7_new_second+f9_new;
a=ft/m;
figure()
plot(a)
title('Compression')
%  plot(f3_new_first)



t1=0.9373*10^(-6);
T=20*10^(-6);
m=0.62;
Ap=12.5*(-(1-exp(-1333.333*(T-t1)))/(1-exp(-1333.33*t1-1333.33*(T-t1))));
Cp=12.5*(1+(exp(-1333.33*t1)*((exp(-1333.33*(T-t1))-1)/(1-exp(-1333.33*t1-1333.33*(T-t1))))));
An=12.5*((1-exp(-1333.33*(T-t1)))/(1-exp(-1333.33*t1-1333.33*(T-t1))));
Cn=12.5*((((1-exp(-1333.33*(T-t1)))*exp(-1333.33*t1))/(1-exp(-1333.33*t1-1333.33*(T-t1))))-1);
iav_p=-37.5*(Ap*exp(-1333.33*t1)-(Ap+Cp)+Cp*exp(-1333.33*(T-t1))-16666.625*t1);
iav_n=-37.5*(Ap*exp(-1333.33*t1)-(Ap+Cp)+Cp*exp(-1333.33*(T-t1))+16666.625*t1);
F_fk=1.2152;
% F_fk=1.8228;
% F_fk=2.4304;
% F_fk=3.038;
n=1;
 for t=0:0.001:1
 x(n)=9.25-9.25*cos(4*pi*t);
   if (x(n)<=18.5)&&(x(n)>17.6)
    R3(n)=11.5-(x(n));
    R1(n)=23.7-(x(n));
    g1(n)=inv(diag(diag((R1(n)')*(R1(n)))));
    g3(n)=inv(diag(diag((R3(n)')*(R3(n)))));
    B_total_3(n)=diag(2.4+((0.4*pi*iav_p)*g3(n)));
    B_total_1(n)=diag(2.4-((0.4*pi*iav_n)*g1(n)));
    k3(n)=diag(B_total_3(n))/(diag(B_total_3(n))+diag(B_total_1(n)));
    k1(n)=diag(B_total_1(n))/(diag(B_total_3(n))+diag(B_total_1(n)));
    f5(n)=(k3(n)+k1(n))*(21.3376*iav_n);

   elseif (x(n)<17.6)&&(x(n)>11.5)
    R4(n)=5.4-(x(n));
    R1(n)=23.7-(x(n));
    g1(n)=inv(diag(diag((R1(n)')*(R1(n)))));
    g4(n)=inv(diag(diag((R4(n)')*(R4(n)))));
    B_total_4(n)=diag(2.4+((0.4*pi*iav_p)*g4(n)));
    B_total_1(n)=diag(2.4-((0.4*pi*iav_n)*g1(n)));
    k4(n)=diag(B_total_4(n))/(diag(B_total_4(n))+diag(B_total_1(n)));
    k1(n)=diag(B_total_1(n))/(diag(B_total_4(n))+diag(B_total_1(n)));
    f4(n)=(k4(n)+k1(n))*(21.3376*iav_n);

   elseif (x(n)<11.5)&&(x(n)>9.25)
    R2(n)=17.6-(x(n));
    R1(n)=23.7-(x(n));
    g1(n)=inv(diag(diag((R1(n)')*(R1(n)))));
    g2(n)=inv(diag(diag((R2(n)')*(R2(n)))));
    B_total_2(n)=diag(2.4-((0.4*pi*iav_n)*g2(n)));
    B_total_1(n)=diag(2.4-((0.4*pi*iav_n)*g1(n)));
    k2(n)=diag(B_total_2(n))/(diag(B_total_2(n))+diag(B_total_1(n)));
    k1(n)=diag(B_total_1(n))/(diag(B_total_2(n))+diag(B_total_1(n)));
    f3(n)=(k2(n)+k1(n))*(21.3376*iav_n);

  elseif (x(n)<9.25)&&(x(n)>5.4)
    R2(n)=17.6-(x(n));
    R1(n)=23.7-(x(n));
    g1(n)=inv(diag(diag((R1(n)')*(R1(n)))));
    g2(n)=inv(diag(diag((R2(n)')*(R2(n)))));
    B_total_2(n)=diag(2.4+((0.4*pi*iav_p)*g2(n)));
    B_total_1(n)=diag(2.4+((0.4*pi*iav_p)*g1(n)));
    k2(n)=diag(B_total_2(n))/(diag(B_total_2(n))+diag(B_total_1(n)));
    k1(n)=diag(B_total_1(n))/(diag(B_total_2(n))+diag(B_total_1(n)));
    f2(n)=(k2(n)+k1(n))*(21.3376*iav_p);
    
  elseif (x(n)>=0)&&(x(n)<5.4)
% R1(n)=23.7-(x(n));
    R2(n)=17.6-(x(n));
    R3(n)=11.5-(x(n));
% R4(n)=5.4-(x(n));
    g2(n)=inv(diag(diag((R2(n)')*(R2(n)))));
    g3(n)=inv(diag(diag((R3(n)')*(R3(n)))));
    B_total_2(n)=diag(2.4+((0.4*pi*iav_p)*g2(n)));
    B_total_3(n)=diag(2.4+((0.4*pi*iav_p)*g3(n)));
    k2(n)=diag(B_total_2(n))/(diag(B_total_2(n))+diag(B_total_3(n)));
    k3(n)=diag(B_total_3(n))/(diag(B_total_2(n))+diag(B_total_3(n)));
    f1(n)=(k2(n)+k3(n))*(21.3376*iav_p);
   end
 n=n+1;
end

f1(1,(length(f1))+1:(length(x)))=zeros(1,(length(x))-(length(f1)));
f2(1,(length(f2))+1:(length(x)))=zeros(1,(length(x))-(length(f2)));
f3(1,(length(f3))+1:(length(x)))=zeros(1,(length(x))-(length(f3)));
f4(1,(length(f4))+1:(length(x)))=zeros(1,(length(x))-(length(f4)));
f5(1,(length(f5))+1:(length(x)))=zeros(1,(length(x))-(length(f5)));
ft=2*(f1+f3+f4+f2+f5);
figure()
a=ft;
plot(a)
title('Expansion')