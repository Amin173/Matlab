clc 
clear
syms V1 V2 V3 V4 V5 V6 Vm e2 e4 h1 h2 h3 h4 w g
% h1=2.5e-5;
% h2=11e-5;
% h3=10e-5;
% h4=0.02;
% g=0.4e-3;
% w=0.4e-3;
e1=3.4;
e3=1;
% V=[0 0 0 0 0 0 0;V7 V6 V7 V8 V9 V10 V9; V12 V11 V12 5000 5000 5000 5000; V17 V16 V17 V18 V19 V20 V19;V22 V21 V22 V23 V24 V25 V24;V27 V26 V27 V28 V29 V30 V29;V32 V31 V32 V33 V34 V35 V34;V27 V36 V37 V38 V39 V40 V39; 0 0 0 0 0 0 0]; 
% V=[0 0 0 0 0;5000 V4 5000 5000 5000; V8 V7 V8 V9 V8;V11 V10 V11 V12 V11; V14 V13 V14 V15 V14];
V=[0 0 0 0 0 0 0;5000 5000 5000 Vm 5000 5000 5000;V2 V1 V2 V3 V2 V1 V2;V5 V4 V5 V6 V5 V4 V5;0 0 0 0 0 0 0];
y=[h1;h2;h3;h4];
x=[g/2;w/2;g/2;g/2;w/2;g/2];
e=[e2;e1;e2;e3;e4;e3];
[a,b]=size(V);
E(1:a,1:b)=poly2sym(0);

for i=2:a-1
    for j=2:b-1
    E(i,j)=e(i)*(x(j)+x(j-1))*(V(i+1,j)-V(i,j))/y(i)+e(i-1)*(x(j)+x(j-1))*(V(i-1,j)-V(i,j))/y(i-1)+0.5*(e(i-1)+e(i))*(y(i)+y(i-1))*((V(i,j+1)-V(i,j))/x(j))+0.5*(y(i)+y(i-1))*(e(i-1)+e(i))*((V(i,j-1)-V(i,j))/x(j-1));
    end 
end

Eq(1)=E(2,4)==0;
Eq(2)=E(3,2)==0;
Eq(3)=E(3,3)==0;
Eq(4)=E(3,4)==0;
Eq(5)=E(4,2)==0;
Eq(6)=E(4,3)==0;
Eq(7)=E(4,4)==0;
% sol=solve(Eq,[V1 V2 V3 V4 V5 V6 Vm]);
% simplify(sol.Vm,'Steps', 300)