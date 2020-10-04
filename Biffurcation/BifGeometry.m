% clc
% clear
% close all
% 
% %% Parameter initialization
% B = [0.0001,0.0001,0.0001];
% va=[1 0 0]+B;
% vb=[0 1 0]+B;
% vc=[-1 -1 0.5]+B;
% 
% 
% Ra=0.3;
% Rb=0.3;
% Rc=0.3;
function P = BifGeometry(va,vb,vc,Ra,Rb,Rc,B)

va=va/norm(va);
vb=vb/norm(vb);
vc=vc/norm(vc);
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

S=S+B;
n= cross((S(3,:)-S(1,:))',(S(2,:)-S(1,:))');
A=mean(S,2);
cp1= B' + n*(Ra+Rb+Rc)/3;
cp2= B' - n*(Ra+Rb+Rc)/3;
%% Vertex coordinates
p1=S(1,:);
p2=S(2,:);
p3=S(3,:);

SLocal=[p1;p2;p3];
CLocal = [cp1,cp2];

P.Spts = SLocal;
P.Cpts = CLocal;

end

% % %%
% x=linspace(B(1),va(1),10000);
% y=linspace(B(2),va(2),10000);
% z=linspace(B(3),va(3),10000);
% r = Ra;
% [X,Y,Z] = tubeplot(x,y,z,r,1,100);
% surf(X,Y,Z,'EdgeColor','none','FaceAlpha',0.01,'FaceColor','b')
% %%
% x=linspace(B(1),vb(1),1000);
% y=linspace(B(2),vb(2),1000);
% z=linspace(B(3),vb(3),1000);
% r = Rb;
% [X,Y,Z] = tubeplot(x,y,z,r,1,100);
% surf(X,Y,Z,'EdgeColor','none','FaceAlpha',0.01,'FaceColor','b')
% %%
% x=linspace(B(1),vc(1),1000);
% y=linspace(B(2),vc(2),1000);
% z=linspace(B(3),vc(3),1000);
% r = Rc;
% [X,Y,Z] = tubeplot(x,y,z,r,10,100);
% surf(X,Y,Z,'EdgeColor','none','FaceAlpha',0.01,'FaceColor','b')