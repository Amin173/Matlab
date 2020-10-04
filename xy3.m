function d=xy3(q)
global l1 l3 l2
x1=l1*cos(q(1));
x2=x1+l2*cos(q(1)+q(2));
x3=x2+l3*cos(q(1)+q(2)+q(3));
y1=l1*sin(q(1));
y2=y1+l2*sin(q(1)+q(2));
y3=y2+l3*sin(q(1)+q(2)+q(3));
d=[x1 x2 x3;y1 y2 y3];
end