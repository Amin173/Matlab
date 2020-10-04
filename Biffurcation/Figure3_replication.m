clc
clear
close all

%% Parameter initialization
B=[0 0 0];
va=[-2 -2 -2];
vb=[8 9 7];
vc=[4 0 1];
va=va/norm(va);
vb=vb/norm(vb);
vc=vc/norm(vc);

Ra=0.3;
Rb=0.3;
Rc=0.3;
%% Computing alphas
alpha_ab=atan2d(norm(cross(va,vb)),dot(va,vb));
alpha_bc=atan2d(norm(cross(vb,vc)),dot(vb,vc));
alpha_ac=atan2d(norm(cross(va,vc)),dot(va,vc));
alpha_vector=[alpha_ab;alpha_bc;alpha_ac];
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
for i=1:length(alpha_vector)
    if (alpha_vector(i)<=90)||((360-alpha_vector(i))<=90)
S(i,:)=k(i,:)*(R(i,1)/sin(atan2(R(i,1),R(i,2))));
    else
        S(i,:)=k(i,:)*(sum(R(i,:))/2);
    end
end
n= cross((S(3,:)-S(1,:))',(S(2,:)-S(1,:))');
S=S+B;
A=mean(S,2);
cp1= B' + n*(Ra+Rb+Rc)/3;
cp2= B' - n*(Ra+Rb+Rc)/3;
% cp1 = B' + n*(Ra+Rb+Rc)/3;
% cp2 = B' - n*(Ra+Rb+Rc)/3;
%% Vertex coordinates
p1=S(1,:);
p2=S(2,:);
p3=S(3,:);
%% Plots
figure('color','w')
h=patch('Faces',1:3,'Vertices',[p1;p2;p3]);
set(h,'FaceColor','b','EdgeColor','k','LineWidth',1,'FaceAlpha',0.3)
axis equal vis3d
view([60 60])
xlabel('x','FontSize',10)
ylabel('y','FontSize',10)
zlabel('z','FontSize',10)

hold on 
plot3([B(1),va(1)+B(1)],[B(2),va(2)+B(2)],[B(3),va(3)+B(3)],':','LineWidth',2,'Color','r')
plot3([B(1),vb(1)+B(1)],[B(2),vb(2)+B(2)],[B(3),vb(3)+B(3)],':','LineWidth',2,'Color','r')
plot3([B(1),vc(1)+B(1)],[B(2),vc(2)+B(2)],[B(3),vc(3)+B(3)],':','LineWidth',2,'Color','r')

% cylinder2(Ra,va,10)
% cylinder2(Rb,vb,10)
% cylinder2(Rc,vc,10)

plot3([cp1(1) cp2(1)],[cp1(2),cp2(2)],[cp1(3),cp2(3)],'LineWidth',4,'LineStyle',':')
scatter3(cp1(1),cp1(2),cp1(3),40,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75])
    scatter3(cp2(1),cp2(2),cp2(3),40,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75])
        scatter3(B(1),B(2),B(3),40,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k')
view(-30,10)
text(cp1(1)+0.1,cp1(2),cp1(3),'CP_1','FontWeight','bold','FontSize',10)
text(cp2(1)+0.1,cp2(2),cp2(3),'CP_2','FontWeight','bold','FontSize',10)
text(B(1)+0.1,B(2),B(3),'B','FontWeight','bold','FontSize',10)
% % %%
x=linspace(B(1),B(1)+va(1),1000);
y=linspace(B(2),B(2)+va(2),1000);
z=linspace(B(3),B(3)+va(3),1000);
r = Ra;
[X,Y,Z] = tubeplot(x,y,z,r,1,100);
surf(X,Y,Z,'EdgeColor','none','FaceAlpha',0.01,'FaceColor','b')
%%
x=linspace(B(1),vb(1)+B(1),1000);
y=linspace(B(2),vb(2)+B(2),1000);
z=linspace(B(3),vb(3)+B(3),1000);
r = Rb;
[X,Y,Z] = tubeplot(x,y,z,r,1,100);
surf(X,Y,Z,'EdgeColor','none','FaceAlpha',0.01,'FaceColor','b')
%%
x=linspace(B(1),vc(1)+B(1),1000);
y=linspace(B(2),vc(2)+B(2),1000);
z=linspace(B(3),vc(3)+B(3),1000);
r = Rc;
[X,Y,Z] = tubeplot(x,y,z,r,10,100);
surf(X,Y,Z,'EdgeColor','none','FaceAlpha',0.01,'FaceColor','b')