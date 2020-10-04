%% Problem 1
clc
clear
syms x1 x2
z=[atan(x1/x2);atan((x1-20)/x2)];
simplify(diff(z,x1))
simplify(diff(z,x2))
%%
clc
clear
th1=46.2*pi/180 ;
th2= 1.4*pi/180;

x1=20;
x2=20;

x_bar=[x1;x2];
Z_bar=[th1;th2];

x_star=x_bar;
H=[x2/(x1^2 + x2^2)  -x1/(x1^2 + x2^2);
   1/(x2*((x1 - 20)^2/x2^2 + 1))   -(x1 - 20)/(x1^2 - 40*x1 + x2^2 + 400)];
h_star=[atan(x1/x2);atan((x1-20)/x2)];
V=eye(2)*(0.3*pi/180)^2;
Z_star=Z_bar-h_star+H*x_star;
x_hat=((H'*inv(V)*H)\H')*(V\Z_star);
P_hat=inv(H'*inv(V)*H);
eps=1e-12;
while((x_hat-x_star)'*(P_hat^-1)*(x_hat-x_star)>eps)
x_star=x_hat;
x1=x_star(1);
x2=x_star(2);
H=[x2/(x1^2 + x2^2)  -x1/(x1^2 + x2^2);
   1/(x2*((x1 - 20)^2/x2^2 + 1))   -(x1 - 20)/(x1^2 - 40*x1 + x2^2 + 400)];
h_star=[atan(x1/x2);atan((x1-20)/x2)];

Z_star=Z_bar-h_star+H*x_star;
x_hat=((H'*inv(V)*H)\H')*(V\Z_star);
P_hat=inv(H'*inv(V)*H);
disp((x_hat-x_star)'*P_hat^-1*(x_hat-x_star));
end
[V,D]=eig(P_hat)
%% Problem 2
clc
clear
dt_MLS=1;
dt_DME=70e-6;
V=[0.01^2 0;0 (0.3e-6)^2];
w=19;
z=[dt_MLS;dt_DME]-[40/w+0.1;5e-6];


c=3e5;
H=[-2/w 0;0 2/c];
x_hat=(H'*(V\H))\(H'*(V\z))
P_hat=inv(H'*(V\H))

%% estimating the height of the aircraft
th=x_hat(1);
r=x_hat(2);
H=[r*cosd(th) sind(th);-r*sind(th) cosd(th)];
P_hat=inv(H'*inv(V)*H)

%% Problem 3
clc
clear
syms x1 x2
h=[sqrt(x1^2+x2^2);sqrt((x1-1000)^2+x2^2);sqrt(x1^2+(x2-1000)^2)];
simplify(diff(h,x1))
simplify(diff(h,x2))
%%
clc
clear
sigma_v=0.1;
V=eye(3)*(sigma_v)^2;
 x1=-550:20:1550;
 x2=-550:20:1550;
 T=zeros(length(x1));
for i=1:length(x1)
for j=1:length(x2)

H=[                x1(i)/(x1(i)^2 + x2(i)^2)^(1/2)   x2(i)/(x1(i)^2 + x2(i)^2)^(1/2);
 (x1(i) - 1000)/((x1(i) - 1000)^2 + x2(i)^2)^(1/2)   x2(i)/((x1(i) - 1000)^2 + x2(i)^2)^(1/2);
          x1(i)/((x2(i) - 1000)^2 + x1(i)^2)^(1/2)  (x2(i) - 1000)/((x2(i) - 1000)^2 + x1(i)^2)^(1/2)];
      P_hat=inv(H'/V*H);
      T(i,j)=trace(P_hat);
end
end
contour(x1',x2',T)
title('Covariance trace of the user’s position')
zlabel('Covariance trace')
ylabel('y [km]')
xlabel('x [km]')

%% Problem 4

clc
clear

ya=10;
yb=7;
xb=50;
l=51;
sigma_v=0.1;
AB=sqrt((ya-yb)^2+xb^2);

th_star=160*pi/180;

h=(l^2-AB^2)/(2*(l+AB*cos(th_star)));
H=AB*(l^2-AB^2)*sin(th_star)/(2*(l+AB*cos(th_star))^2);
V=(sigma_v)^2;
P_hat=inv(H'*inv(V)*H);


% while((th_hat-th_star)'*(P_hat^-1)*(th_hat-th_star)>0.01)
% th_star=th_hat;
% 
% h=(l^2-AB^2)/(2*(l+AB*cos(th_star)));
% H=AB*(l^2-AB^2)*sin(th_star)/(2*(l+AB*cos(th_star))^2);
% 
% Z_star=Z_bar-h+H*th_star;
% th_hat=((H'*inv(V)*H)\H')*(V\Z_star);
% P_hat=inv(H'*inv(V)*H);
% disp((th_hat-x_star)'*P_hat^-1*(th_hat-x_star));
% end

%%
clc
clear
syms th r
th_p=3.5*pi/180;
h=-r*[cos(th-th_p);sin(th-th_p)];
diff(h,th)
diff(h,r)
%%
clc
clear
l=51;
ya=10;
yb=7;
xb=50;
AB=sqrt((ya-yb)^2+xb^2);

th_star=160*pi/180;
r_star=(l^2-AB^2)/(2*(l+AB*cos(th_star)));
sigma_R=0.1;
sigma_th=sqrt(3.84e-6);
V=[sigma_th^2 0;0 sigma_R^2];


H=[  r_star*sin(th_star - (7*pi)/360)  -cos(th_star - (7*pi)/360)
    -r_star*cos(th_star - (7*pi)/360)  -sin(th_star - (7*pi)/360)];

% h=-r_star*[cos(th_star- (7*pi)/360);sin(th_star- (7*pi)/360)];


P_hat=inv(H'*inv(V)*H);

 