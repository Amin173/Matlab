% This program computes and plots the bezier curve approximation of a set
% of given curves in three dimentions
clc
clear
close all
%% Adding tubeplot function path (tubeplot is used for priliminary visaliztion of the vessels)
addpath('C:\Users\amink\Documents\GitHub\Matlab-files\Tubeplot');
%% Arbtrary input data
%n: number of data points in each segment (for simplicity, n is
%
% considered identical for each curve)
n=100;
[x,y,z]=genData(5,n);
x1 = x(:,1);
x2 = x(:,2);
x3 = x(:,5);
x4 = x(:,4);
x5 = x(:,3);

y1 = y(:,1);
y2 = y(:,2);
y3 = y(:,5);
y4 = y(:,4);
y5 = y(:,3);

z1 = z(:,1);
z2 = z(:,2);
z3 = z(:,5);
z4 = z(:,4);
z5 = z(:,3);
%% Bifurcation constriants enforced
[x3(1),y3(1),z3(1)]=deal(x1(end),y1(end),z1(end));
[x2(1),y2(1),z2(1)]=deal(x1(end),y1(end),z1(end));
[x4(1),y4(1),z4(1)]=deal(x2(end),y2(end),z2(end));
[x5(1),y5(1),z5(1)]=deal(x2(end),y2(end),z2(end));
%% Defining the data set
data(:,:,1) = [x1,y1,z1];
data(:,:,2) = [x2,y2,z2];
data(:,:,3) = [x3,y3,z3];
data(:,:,4) = [x4,y4,z4];
data(:,:,5) = [x5,y5,z5];
%% Compute the number of segments
bezierPieces = size(data,3);
%% Compute endpoints
Ra = 0.25*[1 1];
Rb = 0.25*[1 1];
Rc = 0.25*[1 1];

for i=1:size(data,3)
    P(2*i-1,:)=data(1,:,i);
    P(2*i,:)=data(end,:,i);
    %     hold on
    %     figure(1)
    % plot3(data(:,1,i),data(:,2,i),data(:,3,i),'r:','LineWidth',2)
    %     [X,Y,Z] = tubeplot(data(:,1,i),data(:,2,i),data(:,3,i),Ra(1),1,100);
    %     s=surf(X,Y,Z,'EdgeColor','k','EdgeAlpha',0.01,'LineWidth',0.001,'FaceAlpha',0.001);
    %     s.FaceColor='r';
end

Px=zeros(bezierPieces,4);
Py=zeros(bezierPieces,4);
Pz=zeros(bezierPieces,4);
for i = 1:bezierPieces
    nSplinePoints = size(data,1);
    t=setTValueForPts(data(:,:,i),nSplinePoints);
    A0 = zeros(length(t),4);
    A=zeros(length(t),2);
    
    for j = 1:length(t)
        
        b0(j) = (1-t(j)).^3;
        b1(j) = 3*t(j).*(1-t(j)).^2;
        b2(j) = 3*(1-t(j)).*t(j).^2;
        b3(j) = t(j).^3;
        
        A(j,:) = [b1(j) b2(j)];
        A0(j,:) = [b0(j) b1(j) b2(j) b3(j)];
    end
    
    
    xPos = data(:,1,i);
    yPos = data(:,2,i);
    zPos = data(:,3,i);
    
    bVectSingleBezierX = xPos-b0'*P(2*i-1,1)-b3'*P(2*i,1);
    bVectSingleBezierY = yPos-b0'*P(2*i-1,2)-b3'*P(2*i,2);
    bVectSingleBezierZ = zPos-b0'*P(2*i-1,3)-b3'*P(2*i,3);
    
    Cx = A\bVectSingleBezierX;
    Cy = A\bVectSingleBezierY;
    Cz = A\bVectSingleBezierZ;
    
    Pvaluex1 = [P(2*i-1,1) Cx(1) Cx(2) P(2*i,1)];
    Pvaluey1 = [P(2*i-1,2) Cy(1) Cy(2) P(2*i,2)];
    Pvaluez1 = [P(2*i-1,3) Cz(1) Cz(2) P(2*i,3)];
    
    %save control points and endpoints
    Px(i,:)=Pvaluex1;
    Py(i,:)=Pvaluey1;
    Pz(i,:)=Pvaluez1;
    A0s{i}=A0;
    
    finalx1 = A0*Pvaluex1';
    finaly1 = A0*Pvaluey1';
    finalz1 = A0*Pvaluez1';
    
    % plot
    figure(1)
    hold on
    plot3(finalx1,finaly1,finalz1,'LineWidth',2);
    
    % scatter3(Cx(1),Cy(1),Cz(1),'k');  scatter3(Cx(2),Cy(2),Cz(2),'b');
    % text(Cx(1),Cy(1),Cz(1),['C1-',num2str(i)])
    % text(Cx(2),Cy(2),Cz(2),['C2-',num2str(i)])
end
xlabel('x')
ylabel('y')
zlabel('z')
grid on

bifpts = [x1(end),y1(end),z1(end); x2(end),y2(end),z2(end)];
bifspline = [1,2,3;2,4,5];
for i=1:size(bifspline,1)
    C1=[Px(bifspline(i,1),3);Py(bifspline(i,1),3);Pz(bifspline(i,1),3) ];
    P1=[Px(bifspline(i,1),4);Py(bifspline(i,1),4);Pz(bifspline(i,1),4) ];
    va(:,i)=3*C1-3*P1;
end

for i=1:size(bifspline,1)
    C1=[Px(bifspline(i,2),2);Py(bifspline(i,2),2);Pz(bifspline(i,2),2) ];
    P1=[Px(bifspline(i,1),1);Py(bifspline(i,1),1);Pz(bifspline(i,1),1) ];
    vb(:,i)=3*C1-3*P1;
end

for i=1:size(bifspline,1)
    C1=[Px(bifspline(i,3),2);Py(bifspline(i,3),2);Pz(bifspline(i,3),2) ];
    P1=[Px(bifspline(i,3),1);Py(bifspline(i,3),1);Pz(bifspline(i,3),1) ];
    vc(:,i)=3*C1-3*P1;
end
radius = zeros(3,2);
for i = 1:size(bifpts,1)
    B=bifpts(i,:);
    va1 = va(:,i)';
    vb1 = vb(:,i)';
    vc1 = vc(:,i)';
    Ra1 = Ra(i);
    Rb1 = Rb(i);
    Rc1 = Rc(i);
    
    p = BifGeometry(va1,vb1,vc1,Ra1,Rb1,Rc1,B);
    ctrlpts = p.Cpts;
    Spts = p.Spts;
    
    p1=Spts(1,:);
    p2=Spts(2,:);
    p3=Spts(3,:);
    cp1 = ctrlpts(:,1);
    cp2 = ctrlpts(:,2);
    
    f1=BifCurve(cp1',p1,B);
    f2=BifCurve(cp1',p2,B);
    f3=BifCurve(cp1',p3,B);
    f4=BifCurve(cp2',p1,B);
    f5=BifCurve(cp2',p2,B);
    f6=BifCurve(cp2',p3,B);
    
    % plot3(f1(1,:),f1(2,:),f1(3,:),'o')
    % plot3(f2(1,:),f2(2,:),f2(3,:),'o')
    % plot3(f3(1,:),f3(2,:),f3(3,:),'o')
    % plot3(f4(1,:),f4(2,:),f4(3,:),'o')
    % plot3(f5(1,:),f5(2,:),f5(3,:),'o')
    % plot3(f6(1,:),f6(2,:),f6(3,:),'o')
    
    % Plots
    hold on
    h=patch('Faces',1:3,'Vertices',[p1;p2;p3]);
    set(h,'FaceColor','b','EdgeColor','k','LineWidth',1,'FaceAlpha',0.3)
    axis equal vis3d
    xlabel('x','FontSize',10)
    ylabel('y','FontSize',10)
    zlabel('z','FontSize',10)
    
    hold on
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
    %     text(cp1(1)+0.1,cp1(2)-0.1,cp1(3),'CP_1','FontWeight','bold','FontSize',10)
    %     text(cp2(1)+0.1,cp2(2)-0.1,cp2(3),'CP_2','FontWeight','bold','FontSize',10)
    %     text(B(1)+0.1,B(2),B(3),'B','FontWeight','bold','FontSize',10)
    %% Mahsa's function
    controlPoints=ctrlpts';
    separationPoints=Spts;
    bifCoord=B;
    nWalls=100;
    r1=sum(sqrt((controlPoints-B).^2),2);
    r1=r1(1);
    r2=sum(sqrt((separationPoints-B).^2),2);
    radius=[r2, r1*ones(3,1)];
    
    wallDistribution = computeWallDistribution(controlPoints,separationPoints,bifCoord,nWalls);
    
    %------------------- find bifurcation site points -------------------------
    %%following equation have [t 1-t] in their formulae, so type 1 is used
    bifSitPoints1=computeBifSitPointsType1(controlPoints(1,:),separationPoints(1,:),bifCoord,wallDistribution(1,1));
    bifSitPoints3=computeBifSitPointsType1(controlPoints(2,:),separationPoints(2,:),bifCoord,wallDistribution(3,1));
    bifSitPoints8=computeBifSitPointsType1(separationPoints(2,:),controlPoints(2,:),bifCoord,wallDistribution(3,1));
    bifSitPoints5=computeBifSitPointsType1(controlPoints(2,:),separationPoints(3,:),bifCoord,wallDistribution(3,3));
    
    %%following equation have [t 1-t] in their formulae, so type 2 is used
    bifSitPoints2= computeBifSitPointsType1(separationPoints(1,:),controlPoints(2,:),bifCoord,wallDistribution(2,1));
    bifSitPoints4= computeBifSitPointsType1(separationPoints(2,:),controlPoints(1,:),bifCoord,wallDistribution(4,1));
    bifSitPoints7= computeBifSitPointsType1(controlPoints(1,:),separationPoints(2,:),bifCoord,wallDistribution(4,1));
    bifSitPoints6= computeBifSitPointsType1(separationPoints(3,:),controlPoints(1,:),bifCoord,wallDistribution(4,3));
    
    %%arrange bifurcation site points
    bifurSitPoints{1}=[bifSitPoints1;bifSitPoints2;bifSitPoints3;bifSitPoints4];
    bifurSitPoints{2}=[bifSitPoints1;bifSitPoints2;bifSitPoints5;bifSitPoints6];
    bifurSitPoints{3}=[bifSitPoints7;bifSitPoints8;bifSitPoints5;bifSitPoints6];
    %--------------------------------------------------------------------------
    bifSitPoints1 = [bifSitPoints1;separationPoints(1,:)];
    bifSitPoints3 = [bifSitPoints3;separationPoints(2,:)];
    bifSitPoints8 = [bifSitPoints8;controlPoints(2,:)];
    bifSitPoints5 = [bifSitPoints5;separationPoints(3,:)];
    bifSitPoints2 = [bifSitPoints2;controlPoints(2,:)];
    bifSitPoints4 = [bifSitPoints4;controlPoints(1,:)];
    bifSitPoints7 = [bifSitPoints7;separationPoints(2,:)];
    bifSitPoints6 = [bifSitPoints6;controlPoints(1,:)];
    hold on
    scatter3(separationPoints(:,1),separationPoints(:,2),separationPoints(:,3),'filled','b','MarkerEdgeColor','k'); ...
        scatter3(controlPoints(:,1),controlPoints(:,2),controlPoints(:,3),'filled','r','MarkerEdgeColor','k');   ...
        scatter3(bifCoord(:,1),bifCoord(:,2),bifCoord(:,3),'filled','k','MarkerEdgeColor','k');
    %     scatter3(bifSitPoints1(:,1),bifSitPoints1(:,2),bifSitPoints1(:,3))
    %     scatter3(bifSitPoints2(:,1),bifSitPoints2(:,2),bifSitPoints2(:,3))
    %     scatter3(bifSitPoints3(:,1),bifSitPoints3(:,2),bifSitPoints3(:,3))
    %     scatter3(bifSitPoints4(:,1),bifSitPoints4(:,2),bifSitPoints4(:,3))
    %     scatter3(bifSitPoints5(:,1),bifSitPoints5(:,2),bifSitPoints5(:,3))
    %     scatter3(bifSitPoints6(:,1),bifSitPoints6(:,2),bifSitPoints6(:,3))
    %% Saving variables
    p1Sav(i,:)=p1;
    p2Sav(i,:)=p2;
    p3Sav(i,:)=p3;
    CtrlP1(i,:)=cp1';
    CtrlP2(i,:)=cp2';
end

%% Bifurcation point local cooordinate system
%First Axis
[ex,ey,ez] = bifaxis(CtrlP1(1,:),p3Sav(1,:),bifpts(1,:));
%%
hold on
plot3([bifpts(1,1) ex(1)],[bifpts(1,2) ex(2)],[bifpts(1,3) ex(3)],'k-','LineWidth',2);
plot3([bifpts(1,1) ey(1)],[bifpts(1,2) ey(2)],[bifpts(1,3) ey(3)],'k-','LineWidth',2);
plot3([bifpts(1,1) ez(1)],[bifpts(1,2) ez(2)],[bifpts(1,3) ez(3)],'k-','LineWidth',2);

text(ex(1)+0.1,ex(2),ex(3),'ex','FontSize',12);
text(ey(1)+0.1,ey(2),ey(3),'ey','FontSize',12);
text(ez(1)+0.1,ez(2),ez(3),'ez','FontSize',12);
text(p3Sav(1,1)+0.1,p3Sav(1,2)+0.1,p3Sav(1,3),'SP3-Bif_1','FontSize',12,'Color','r');
%%
splineCtrlPts=[Px(2,:);Py(2,:);Pz(2,:)];
crossSection=Circle3D(Ra(2),splineCtrlPts(:,1),ez-bifpts(1,:));
NSegments=500;
q=SurfMesh(0,0.5,NSegments,ex-bifpts(1,:),ey-bifpts(1,:),ez-bifpts(1,:),crossSection,splineCtrlPts);
mesh=q.SMesh;
CSNormals=q.Nvector;
for i=1:NSegments
A=mesh{i};
hold on
plot3(A(1,:),A(2,:),A(3,:),'r-');
end

% Second axis
[ex,ey,ez] = bifaxis(CtrlP1(2,:),p1Sav(2,:),bifpts(2,:));
%%
hold on
plot3([bifpts(2,1) ex(1)],[bifpts(2,2) ex(2)],[bifpts(2,3) ex(3)],'k-','LineWidth',2);
plot3([bifpts(2,1) ey(1)],[bifpts(2,2) ey(2)],[bifpts(2,3) ey(3)],'k-','LineWidth',2);
plot3([bifpts(2,1) ez(1)],[bifpts(2,2) ez(2)],[bifpts(2,3) ez(3)],'k-','LineWidth',2);

text(ex(1)+0.1,ex(2),ex(3),'ex','FontSize',12);
text(ey(1)+0.1,ey(2),ey(3),'ey','FontSize',12);
text(ez(1)+0.1,ez(2),ez(3),'ez','FontSize',12);
text(p1Sav(2,1)+0.1,p1Sav(2,2)+0.1,p1Sav(2,3),'SP1-Bif_2','FontSize',12,'Color','r');
%%
splineCtrlPts2=flip(splineCtrlPts')';
crossSection=Circle3D(Ra(2),splineCtrlPts2(:,1),ey-bifpts(2,:));
q2=SurfMesh(0,0.5,NSegments,ez-bifpts(2,:),ex-bifpts(2,:),ey-bifpts(2,:),crossSection,splineCtrlPts2);
mesh2=q2.SMesh;
CSNormals2=q2.Nvector;
for i=1:NSegments
A=mesh2{i};
hold on
plot3(A(1,:),A(2,:),A(3,:),'r-');
end
%%
firstCSection=mesh{end};
secondCSection=mesh2{end};
[r1,center,r2,ii,jj] = closestPts(firstCSection,secondCSection);

%% Computing the minimum roll angle
minRoll=computeAngle(r1,center,r2);
%% ploting the two midle cross sections from t=0-0.5 and t=1-0.5
figure(2)
hold on   
plot3(firstCSection(1,:),firstCSection(2,:),firstCSection(3,:),'r-');
plot3(secondCSection(1,:),secondCSection(2,:),secondCSection(3,:),'b-');
plot3([center(1) r1(1)],[center(2) r1(2)],[center(3) r1(3)]);
plot3([center(1) r2(1)],[center(2) r2(2)],[center(3) r2(3)]);
%% rotating the cross sections

for i=1:NSegments
A=mesh{i};
B=mesh2{i};
Va=CSNormals{i};
Vb=CSNormals2{i};
A=fixRoll(A,minRoll/2,Va,ii,jj);
B=fixRoll(B,-minRoll/2,Vb,jj,ii);
mesh{i}=A;
mesh2{i}=B;
end
%%
firstCSection=mesh{end};
secondCSection=mesh2{end};
r1=firstCSection(:,ii);
r2=secondCSection(:,jj);
minRoll2=computeAngle(r1,center,r2); % for checking the roll elimination sucsess

figure(3)
hold on
plot3(firstCSection(1,:),firstCSection(2,:),firstCSection(3,:),'r','LineWidth',2);
plot3(secondCSection(1,:),secondCSection(2,:),secondCSection(3,:),'b','LineWidth',2);
plot3([center(1) r1(1)],[center(2) r1(2)],[center(3) r1(3)]);
plot3([center(1) r2(1)],[center(2) r2(2)],[center(3) r2(3)]);
%% Ploting new cross sections
figure(4)
hold on
for i=1:NSegments
A=mesh{i};
hold on
plot3(A(1,:),A(2,:),A(3,:),'r-');
end
for i=1:NSegments
A=mesh2{i};
hold on
plot3(A(1,:),A(2,:),A(3,:),'r-');
end
