clc
clear
l1=1;l2=1;
lc1=0.*l1;
lc2=0.5*l2;
% J=[-l1*sin(x(1))-lc2*sin(x(1)+x(2)) -lc2*sin(x(1)+x(2));l1*cos(x(1))+lc2*cos(x(1)+x(2)) lc2*cos(x(1)+x(2))];
y1=0.2+0.25;
y2=0.2-0.25:0.05:0.2+0.25;
y3=0.2-0.25;
x1=0.5-0.25;
x2=0.5-0.25:0.05:0.5+0.25;
x3=0.5+0.25;
k=0;
for j=1:4
switch j
    case 1
        y=y1;
        x=x2;
    case 2
        y=y2;
        x=x3;
    case 3
        y=y3;
        x=x2;
    case 4
        y=y2;
        x=x1;
end
a=max([length(x),length(y)]);
if length(x)<a
    x=x*y./y;
else
    y=y*x./x;
end
for i=1:a
    k=k+1;
c2(k)=((x(i).^2+y(i).^2-l1^2-l2^2)/(2*l1*l2));
s2(k)=(1-c2(k).^2).^0.5;
D=inv([l1+l2*c2(k) -l2*s2(k);l2*s2(k) l1+l2*c2(k)])*[x(i);y(i)];
c1(k)=D(1);
s1(k)=D(2);
theta1(k)=atan2(c1(k),s1(k));
theta2(k)=atan2(c2(k),s2(k));
hold on
plot(theta1(k),theta2(k),'.','linewidth',3)
end
end
% plot(theta1,theta2);
axis equal
xlim([1 3.5])
ylim([-2 0])
xlabel('theta1(rad)')
ylabel('theta2(rad)')