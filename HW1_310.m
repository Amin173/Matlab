clc
clear
close all
%% HW one of BioE 310
A = [1, 2, 2^2, 2^3;
    1, 3, 3^2, 3^3;
    1, 4, 4^2, 4^3;
    1, 6, 6^2, 6^3];
b = [0; 23; 3; 40];
coef = A\b;
coef = flipud(coef);
x1 = [2,3,4,6];
y = [0,23,3,40];
x = 2:0.1:6;
y1 = polyval(coef,x);
k=polyder(coef);
m=polyval(k,x1);
y0=polyval(coef,x1);
%% Lagrangian Form
% l0 = ((x1-x(2)).*(x1-x(3)).*(x1-x(4)))/((x(1)-x(2)).*(x(1)-x(3)).*(x(1)-x(4)));
% l1 = ((x1-x(1)).*(x1-x(3)).*(x1-x(4)))/((x(2)-x(1)).*(x(2)-x(3)).*(x(2)-x(4)));
% l2 = ((x1-x(1)).*(x1-x(2)).*(x1-x(4)))/((x(3)-x(1)).*(x(3)-x(2)).*(x(3)-x(4)));
% l3 = ((x1-x(1)).*(x1-x(2)).*(x1-x(3)))/((x(4)-x(1)).*(x(4)-x(2)).*(x(4)-x(3)));
% or this way
l0 = ((x.^3)-(13*x.^2)+(54.*x)-72)./(-8);
l0_tan=((3*x.^2)-(26*x)+54)./(-8);
l1 = ((x.^3)-(12*x.^2)+(44.*x)-48)./3;
l1_tan=((3*x.^2)-(24*x)+44)./3;
l2 = ((x.^3)-(11*x.^2)+(36.*x)-36)./(-4);
l2_tan=((3*x.^2)-(22*x)+36)./(-4);
l3 = ((x.^3)-(9*x.^2)+(26.*x)-24)./24;
l3_tan=((3*x.^2)-(18*x)+26)./24;
%% finding extreme points
[xe0,ye0] = find_extreme(x,l0);
[xe1,ye1] = find_extreme(x,l1);
[xe2,ye2] = find_extreme(x,l2);
[xe3,ye3] = find_extreme(x,l3);
%% plotting
figure(1)
plot(x,y1,'*-')
hold on
f = l1.*23+l2.*3+l3.*40;
plot(x,f,'k')
hold on
plot(x1,y,'ro')
dx=0.25;
for i=1:length(x1)
    plot([x1(i)-dx x1(i)+dx],[m(i)*(-dx)+y0(i) m(i)*dx+y0(i)],'c','LineWidth',1.5)
end
figure(2)
hold on
plot(x,l0,'p', x,l1,'b', x,l2,'k',x,l3,'r')
xe=[xe0;xe1;xe2;xe3];
ye=[ye0;ye1;ye2;ye3];
for i=1:size(xe,1)
plot([xe(i)-0.25 xe(i)+0.25],[ye(i) ye(i)],'c','LineWidth',1.5)
end
y3 = [0,0,0,0];
plot(x1,y3,'*')