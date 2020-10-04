clc
clear

r0 = 1.5 * 2.54;
l0 = 6* 2.54;
Rmin = l0/(pi);
R = Rmin+2:1:round(Rmin)+100;

for j = 1:length(R)
theta = pi-l0/R(j);
% P0xu(j) = -(R(j)+r0)*cos(theta/2);
% P2xu(j) = -P0xu(j);
% P0yu(j)= (R(j)+r0)*sin(theta/2);
% P2yu(j) = P0yu(j);
% P1xu(j) = 0;
P0xu(j) = -(R(j)+2*r0)*cos(theta/2);
P2xu(j) = -P0xu(j);
P0yu(j)= (R(j)+2*r0)*sin(theta/2);
P2yu(j) = P0yu(j);
P1xu(j) = 0;
df = @(t, P1y)sqrt((2*(1-t)*(P1xu(j)-P0xu(j)) + 2*t*(P2xu(j)-P1xu(j))).^2 + (2*(1-t)*(P1y-P0yu(j)) + 2*t*(P2yu(j)-P1y)).^2);
n = 100;
dt = 1/n;
t = 0:dt:1;
f = @(P1y) abs(dt * 0.5 * sum(df(t(1:end-1), P1y) +  df(t(2:end), P1y))-l0);
options = optimset('TolX',1e-10);
P1yu(j) = fminbnd(f, R(j), 5*R(j), options);
% end
% 
% for j = 1:length(R)
% theta = pi-l0/R(j);
% P0xl(j) = -(R(j)-r0)*cos(theta/2);
% P2xl(j) = -P0xl(j);
% P0yl(j)= (R(j)-r0)*sin(theta/2);
% P2yl(j) = P0yl(j);
% P1xl(j) = 0;
P0xl(j) = -(R(j))*cos(theta/2);
P2xl(j) = -P0xl(j);
P0yl(j)= (R(j))*sin(theta/2);
P2yl(j) = P0yl(j);
P1xl(j) = 0;
df = @(t, P1y)sqrt((2*(1-t)*(P1xl(j)-P0xl(j)) + 2*t*(P2xl(j)-P1xl(j))).^2 + (2*(1-t)*(P1y-P0yl(j)) + 2*t*(P2yl(j)-P1y)).^2);
n = 1000;
dt = 1/n;
t = 0:dt:1;
f = @(P1y) abs(dt * 0.5 * sum(df(t(1:end-1), P1y) +  df(t(2:end), P1y))-l0);
options = optimset('TolX',1e-10);
P1yl(j) = fminbnd(f, 2 * R(j) - 0.5 * (P0yl(j)+P2yl(j)), 5*R(j), options);
% P1yl(j) = 2 * R(j) - 0.5 * (P0yl(j)+P2yl(j));
end
%%
t= linspace(0,1,100);
A=zeros(length(t), length(R));
for j= 1:length(R)
Xu = (1-t).^2 * P0xu(j) + 2*(1-t).*t * P1xu(j) + t.^2 * P2xu(j);
Yu = (1-t).^2 * P0yu(j) + 2*(1-t).*t * P1yu(j) + t.^2 * P2yu(j);
Xl = (1-t).^2 * P0xl(j) + 2*(1-t).*t * P1xl(j) + t.^2 * P2xl(j);
Yl = (1-t).^2 * P0yl(j) + 2*(1-t).*t * P1yl(j) + t.^2 * P2yl(j);
A(:,j) = ((Xu-Xl).^2+(Yu-Yl).^2).^0.5;
end

A = A/2;        % converting diameters to radius
p = 2*pi*r0;    % constant perimeter of each cross section
B = abs(2*r0^2 - A.^2).^0.5;    %calculating minor axis of the cross section
CS = pi * A.*B;     %calculating the area

for i = 1:size(CS, 2)
    V(i) = l0*trapz(t, CS(:,i));
end

figure(1)
hold on
plot(R, V,'LineWidth',2)
plot(R, pi*r0^2*l0*ones(length(R),1),'LineWidth',2)
xlabel('Object radius (R) [cm]')
ylabel('Membrane volume [cm^3]')
legend('Deformed membrane volume', 'Undeformed membrane volume')
ylim([min(V)-5, max(V)+10])
%%
clc
for j=1:5:length(R)
Xu = (1-t).^2 * P0xu(j) + 2*(1-t).*t * P1xu(j) + t.^2 * P2xu(j);
Yu = (1-t).^2 * P0yu(j) + 2*(1-t).*t * P1yu(j) + t.^2 * P2yu(j);
Xl = (1-t).^2 * P0xl(j) + 2*(1-t).*t * P1xl(j) + t.^2 * P2xl(j);
Yl = (1-t).^2 * P0yl(j) + 2*(1-t).*t * P1yl(j) + t.^2 * P2yl(j);
X= [];
Y=[];
Z=[];
    for i = 1:length(t)
        pitch = atan2((Yu(i)-Yl(i)),(Xu(i)-Xl(i)))+pi/2;
        [coordinates, h] = ellipse3D(B(i, j), A(i, j), mean([Xu(i),Xl(i)]), mean([Yu(i),Yl(i)]), 0, 100, -pitch, pi/2, 0);
        X = [X;coordinates(1,:)'];
        Y = [Y;coordinates(2,:)'];
        Z = [Z;coordinates(3,:)'];
    end
X=reshape(X,100,100);
Y=reshape(Y,100,100);
Z=reshape(Z,100,100);
figure()
hold on
C = abs(Z);
s = surf(X,Y,Z,C,'FaceAlpha',0.7);
s.EdgeColor = 'none';
% s.CData = C;
% plot(Xl,Yl,'*')
% plot(Xu,Yu,'*')
% plot3(X, Y, Z)
[xt, yt, zt] = cylinder(r0, 1000);
Xt = linspace(-l0/2, l0/2, 2)'.*ones(2, 1001);
Yt = xt + R(j) + r0;
Zt = yt;
s=surf(Xt, Yt, Zt);
s.EdgeColor = 'none';
s.FaceAlpha = 0.3;
s.FaceColor = 'b';
[Xc,Yc,Zc] = cylinder(R(j),100);
s = mesh(Xc,Yc,3*[ones(1,101);-ones(1,101)]);
s.FaceColor = 'k';
zlabel('Z [cm]')
xlabel('X [cm]')
ylabel('Y [cm]')
axis equal
title(['R = ', num2str(round(R(j)))])
end