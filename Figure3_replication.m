clc
clear
close all

%% Parameter initialization
va=[1 0 0];
vb=[-1 1 1];
vc=[-1 -1 10];
va=va/norm(va);
vb=vb/norm(vb);
vc=vc/norm(vc);

Ra=0.4;
Rb=0.5;
Rc=0.6;
%% Computing alphas
alpha_ab=atan2d(norm(cross(va,vb)),dot(va,vb));
alpha_bc=atan2d(norm(cross(vb,vc)),dot(vb,vc));
alpha_ac=atan2d(norm(cross(va,vc)),dot(va,vc));
alpha=[alpha_ab;alpha_bc;alpha_ac];
%% k vectors
kab=Ra*vb+Rb*va;
kab=kab/norm(kab);
kac=Ra*vc+Rc*va;
kac=kac/norm(kac);
kbc=Rb*vc+Rc*vb;
kbc=kbc/norm(kbc);
k=[kab;kac;kbc];
%% 
R=[Ra Rb;Rb Rc;Ra Rc];
for i=1:length(alpha)
    if (alpha(i)<=90)||((360-alpha(i))<=90)
S(i,:)=k(i,:)*(R(i,1)/sin(atan2(R(i,1),R(i,2))));
    else
        S(i,:)=k(i,:)*(sum(R(i,:))/2);
    end
end
n= cross((S(3,:)'-S(1,:)'),(S(2,:)'-S(1,:)'));
B = [0,0,0];
cp1= B + n*(Ra+Rb+Rc)/3;
cp2 = B - n*(Ra+Rb+Rc)/3;

%% Vertex coordinates
p1=S(1,:);
p2=S(2,:);
p3=S(3,:);
%% Plots
figure('color','w')
h=patch('Faces',1:3,'Vertices',[p1;p2;p3]);
set(h,'FaceColor','r','EdgeColor','k','LineWidth',2,'FaceAlpha',0.5)
axis equal vis3d
view([30 30])
xlabel('x','FontSize',20)
ylabel('y','FontSize',20)
zlabel('z','FontSize',20)

hold on 
plot3([B(1),va(1)],[B(2),va(2)],[B(3),va(3)],'LineWidth',2*Ra)
plot3([B(1),vb(1)],[B(2),vb(2)],[B(3),vb(3)],'LineWidth',2*Rb)
plot3([B(1),vc(1)],[B(2),vc(2)],[B(3),vc(3)],'LineWidth',2*Rc)

% cylinder2(Ra,va,10)
% cylinder2(Rb,vb,10)
% cylinder2(Rc,vc,10)

plot3([cp1(1) cp2(1)],[cp1(2),cp2(2)],[cp1(3),cp2(3)],'LineWidth',4,'LineStyle',':')

%%
t=0:0.01:pi;
x=cos(t);
y=sin(t);
z=3*t;
r = sqrt(t)/2;
[X,Y,Z] = tubeplot(x,y,z,r,1,10);
surf(X,Y,Z,'EdgeColor','none')
