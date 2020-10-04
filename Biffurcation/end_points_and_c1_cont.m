clc
clear
close all
%% Bezier function coeficients
t = 0: .001:1;
b0 = (1-t).^3;
b1 = 3*t.*(1-t).^2;
b2 = 3*(1-t).*t.^2;
b3 = t.^3;
%% Control points
P0 = [0,0,0; 1,1,1; 2,2,2; 2,2,2; 1,1,1];
P1 = [0.2,0.8,0.6; 1.2,1.8,1.6; 2.8, 3.5,3; 2.5,1.2,1; 2.5,1.2,0.7];
P2 = [0.8,0.2,0.3; 1.8,1.2,0.9; 3.5,2.2,2.3; 3.5,2.5,2.2; 1.5,1,1.1];
P3 = [1,1,1; 2,2,2; 4,3,1.6; 4,1.5,4; 2,0.5,1.2];
%% Plots
for i=1:size(P0,1)
x = P0(i,1)*b0+P1(i,1)*b1+P2(i,1)*b2+P3(i,1)*b3;
y = P0(i,2)*b0+P1(i,2)*b1+P2(i,2)*b2+P3(i,2)*b3;
z = P0(i,3)*b0+P1(i,3)*b1+P2(i,3)*b2+P3(i,3)*b3;
hold on
figure(1)
plot3(x,y,z)
r=0.025;
[X,Y,Z] = tubeplot(x,y,z,r,1,100);
if i==1
    [X1,Y1,Z1]=deal(X,Y,Z);
else
[X1,Y1,Z1]=deal([X1;X],[Y1;Y],[Z1;Z]);
end
s=surf(X,Y,Z,'EdgeColor','k','EdgeAlpha',0.5,'LineWidth',0.001,'FaceAlpha',1);
s.FaceColor='r';

end
view([40 40])