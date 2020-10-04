clc
clear
global l1 l3 l2

l1=1;
l2=1;
l3=1;
qs=[0;0;0];
qf=(pi/2)*[1;1;1];
obstacle{1}=[1 1.5 1.5 1;2 2 2.5 2.5];
obstacle{2}=[1.5 2 2 1.5;1 1 1.5 1.5]+1;
obstacle{3}=[2.5 3 3 2.5;1 1 1.5 1.5]+1;
eta=50;
zeta=1;
r0=0.2;
q=qs;
qt=qs;
d=4;
for k=1:30
Fatt=-zeta*(xy3(q)-xy3(qf)).*(dist(xy3(qf)-xy3(q))<d)-d*zeta*(xy3(q)-xy3(qf))./dist(xy3(qf)-xy3(q)).*dist(xy3(qf)-xy3(q))>=d;
for i=1:3
    A=[[0;0] xy3(q)];
    P=A(:,i);
    Q=A(:,i+1);
for j=1:3
    O=obstacle{j};
    r=rho(Q,O);
    BB=bmin(O,P,Q);
    b=BB(1,1);
    v=BB(:,2);
   
Frep(:,j)=eta*(1/r-1)^2*(1/r)^2*((Q-v)/norm(Q-v))*(r<r0);
end
    Fr(:,i)=Frep(:,1)+Frep(:,2)+Frep(:,3);
end
J=[-l1*sin(q(1))-l2*sin(q(1)+q(2))-l3*sin(q(1)+q(2)+q(3)), -l2*sin(q(1)+q(2))-l3*sin(q(1)+q(2)+q(3)),-l3*sin(q(1)+q(2)+q(3));
    l1*cos(q(1))+l2*cos(q(1)+q(2))+l3*cos(q(1)+q(2)+q(3)),l2*cos(q(1)+q(2))+l3*cos(q(1)+q(2)+q(3)), l3*cos(q(1)+q(2)+q(3))]; 
F=Fr+Fatt;
tow=J'*F+1e-5;
alpha=0.09;

q=q+alpha*sum(tow,2)/norm(tow);
qt=[qt q];
end

p=xy4(qt);
plot(p(1,:),p(2,:));
xlabel('X(end effector)')
ylabel('Y(end effector)')