%%% MATLAB code for Parametric Meshing of Non-planar Bifurcation
%%% Mahsa 12/10/2013
%%% Linninger research Group


function [ptCoordMx, faceMx, bif, faceMxSize,ptCoordMxSize] = Bifurcation(bzNWK,nwk,nWalls,criteria,wallDensity,fluidDensity)
% clear all;close all;clc
% load ('sample'); axis equal; axis off
%%%%%pre-allocatoin of matrix
bif.circle = [];faceMxSize = []; bif.index =[];
ptCoordMxSize = 0;faceMxSize.in = 0;faceMxSize.inlet = 0;
faceMxSize.sur = 0;faceMxSize.out = 0;
ptCoordMx = zeros(10^7,3);faceMx.in = zeros(10^7,6);
faceMx.inlet = zeros(10^4,6);faceMx.sur = zeros(10^5,6);

for i = 1:(nwk.nBifurcation);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%PRE-PROCESSING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%extracting the bezier cure data for eall three branches of each
    %%bifurcaiton
    bifFaceMx = [];
    spline = bzNWK.obj2msh(i,(2:4));
    bifPoint = nwk.bifNodes(i);
    bifCoord = nwk.ptCoordMx(bifPoint,:);
    d1 = bzNWK.bzSplines(spline(1),1:4);
    d2 = bzNWK.bzSplines(spline(2),1:4);
    d3 = bzNWK.bzSplines(spline(3),1:4);
    Mx=[d1;d2;d3];
    %%%Checking the correct order of spline to make sure all start end of
    %%%bracnh and end with the bifurcation point
    if all(bifCoord == bzNWK.bzPtCrdMx(Mx(1,1),:))==1
        Mx = [fliplr(Mx(1,:));Mx(2:3,:)];
    end
    if all(bifCoord == bzNWK.bzPtCrdMx(Mx(2,4),:))==1
        Mx = [Mx(1,:);fliplr(Mx(2,:));Mx(3,:)];
    end
    if all(bifCoord == bzNWK.bzPtCrdMx(Mx(3,4),:))==1
        Mx = [Mx(1:2,:);fliplr(Mx(3,:))];
    end
    %%%%%extracting C0 and C1 out of bzNWK.bzPtCrdMx
    c0 = [bzNWK.bzPtCrdMx(Mx(1,2),:);bzNWK.bzPtCrdMx(Mx(2,3),:);bzNWK.bzPtCrdMx(Mx(3,3),:)];
    c1 = [bzNWK.bzPtCrdMx(Mx(1,3),:);bzNWK.bzPtCrdMx(Mx(2,2),:);bzNWK.bzPtCrdMx(Mx(3,2),:)];
    pointCoord(1:3,:) = [bzNWK.bzPtCrdMx(Mx(1,1),:);bzNWK.bzPtCrdMx(Mx(2,4),:);bzNWK.bzPtCrdMx(Mx(3,4),:)];
    %     scatter3(c0(:,1),c0(:,2),c0(:,3),'*g');%     scatter3(c1(:,1),c1(:,2),c1(:,3),'*m');
    
    %%%%double-check the order of c0 and c1
    for y = 1:3
        c0Temp = c0(y,:);
        c1Temp = c1(y,:);
        lengthC0 = pdist([c0(y,:);pointCoord(y,:)]);
        lengthC1 = pdist([c1(y,:);pointCoord(y,:)]);
        if lengthC1 < lengthC0
            if y == 1
                c0 = [c1Temp;c0(2,:);c0(3,:)];
                c1 = [c0Temp;c1(2,:);c1(3,:)];
            elseif y == 2
                c0 = [c0(1,:);c1Temp;c0(3,:)];
                c1 = [c1(1,:);c0Temp;c1(3,:)];
            elseif y == 3
                c0 = [c0(1,:);c0(2,:);c1Temp];
                c1 = [c1(1,:);c1(2,:);c0Temp];
            end
        end
    end
    pointTangent=(c0-pointCoord); %velocity at end of the branch
    bifTangent=zeros(3,3);radius = zeros(3,2);
    for u=1:3
        bifTangent(u,:)=-(c1(u,:)-bifCoord);
        tempRadius = bzNWK.diaSpline(spline(u),:);
        tempRadius(tempRadius == 0)=[];
        radius(u,1:2) = tempRadius;
    end
    %%%%bif-bif or terminal-bif detection, calculate length of each branch
    termMx =[];
    for p = 1:size(bzNWK.obj2msh,1)
        if bzNWK.obj2msh(p,1) == 1
            termMx = [termMx;bzNWK.obj2msh(p, 2:3)];
        end
    end
    p = ismember(termMx,spline).*termMx;
    u = p(p ~= 0);
    for t = 1:3
        if spline(1,t)==u;
            splineTerm(1,t)=1;
        else
            splineTerm(1,t)=0;
        end
    end
    arcLen = zeros(3,1);
    for u=1:3
        arcLen(u)=computeArcLength(pointCoord(u,:), c0(u,:),c1(u,:),bifCoord);
        %         scatter3(pointCoord(u,1),pointCoord(u,2),pointCoord(u,3),'filled','b');
    end
    %%%%compute required longitudinal cross-section based on criteria
    %%%%(longitudinal mesh, radius and arc length
    nSteps= computenSteps(arcLen,criteria,radius);
    %%%%optional Draw bezier lines
    for w = 1:3
        centerBezierLine = computeBezier(pointCoord(w,:),c0(w,:),bifCoord,c1(w,:),nSteps(w));%         plot3(centerBezierLine(:,1),centerBezierLine(:,2),centerBezierLine(:,3),':','linewidth',3.3,'color','k');
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Compute seperation points%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [separationPoints ifMoreThan90] = computeSeparationPoints(c1,bifCoord,pointCoord,radius);
%     scatter3(separationPoints(:,1),separationPoints(:,2),separationPoints(:,3),'filled','b');
    ifMoreThan90 = [1 1 1];  %set 0 if you do not want to keep the continuity between branches

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Compute control points%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    controlPoints = computeControlPoints(separationPoints,bifCoord,radius);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%POINTS BETWEEN CONTROL AND SEPERATION POINTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wallDistribution = computeWallDistribution(controlPoints,separationPoints,bifCoord,nWalls);
    
    %------------------- find bifurcation site points -------------------------
    %%following equation have [t 1-t] in their formulae, so type 1 is used
    bifSitPoints1=computeBifSitPointsType1(controlPoints(1,:),separationPoints(1,:),bifCoord,radius,radius(1,:),radius(2,:),wallDistribution(1,1));
    bifSitPoints3=computeBifSitPointsType1(controlPoints(2,:),separationPoints(2,:),bifCoord,radius,radius(1,:),radius(3,:),wallDistribution(3,1));
    bifSitPoints8=computeBifSitPointsType1(separationPoints(2,:),controlPoints(2,:),bifCoord,radius,radius(1,:),radius(3,:),wallDistribution(3,1));
    bifSitPoints5=computeBifSitPointsType1(controlPoints(2,:),separationPoints(3,:),bifCoord,radius,radius(1,:),radius(2,:),wallDistribution(3,3));
    
    %%following equation have [t 1-t] in their formulae, so type 2 is used
    bifSitPoints2= computeBifSitPointsType1(separationPoints(1,:),controlPoints(2,:),bifCoord,radius,radius(1,:),radius(2,:),wallDistribution(2,1));
    bifSitPoints4= computeBifSitPointsType1(separationPoints(2,:),controlPoints(1,:),bifCoord,radius,radius(1,:),radius(3,:),wallDistribution(4,1));
    bifSitPoints7= computeBifSitPointsType1(controlPoints(1,:),separationPoints(2,:),bifCoord,radius,radius(1,:),radius(3,:),wallDistribution(4,1));
    bifSitPoints6= computeBifSitPointsType1(separationPoints(3,:),controlPoints(1,:),bifCoord,radius,radius(1,:),radius(2,:),wallDistribution(4,3));
    
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
    % scatter3(separationPoints(:,1),separationPoints(:,2),separationPoints(:,3),'filled','b','MarkerEdgeColor','k');%     scatter3(controlPoints(:,1),controlPoints(:,2),controlPoints(:,3),'filled','r','MarkerEdgeColor','k') %     scatter3(bifCoord(:,1),bifCoord(:,2),bifCoord(:,3),'filled','k','MarkerEdgeColor','k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Compute wall points on each cross-section%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   previousCir = bif.circle;
    [bif] = computeCirclePoints(pointCoord,separationPoints,pointTangent,radius,nWalls,wallDistribution,controlPoints,bifCoord,bifTangent,spline,previousCir,bif,bifurSitPoints);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Compute all points and INDEXING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [ptCoordMx, faceMx, bif, faceMxSize, ptCoordMxSize] = computePtCoordMx(bif,pointTangent,bifurSitPoints,bifTangent,nWalls,nSteps,c0,c1,pointCoord,bifCoord,ptCoordMx, faceMx,spline,bifFaceMx,splineTerm,ptCoordMxSize,faceMxSize,wallDensity,fluidDensity,ifMoreThan90);
    
end
function [bif] = computeCirclePoints(pointCoord,~,pointTangent,radius,nWalls,~,controlPoints,bifCoord,bifTangent,spline,previousCir,bif,bifurSitPoints)
radius = radius(:,2);
projectedC1_1=computeStartPoint(controlPoints(1,:),bifCoord,bifTangent(1,:),pointCoord(1,:),pointTangent(1,:),radius(1));
projectedC2_1=computeStartPoint(controlPoints(2,:),bifCoord,bifTangent(1,:),pointCoord(1,:),pointTangent(1,:),radius(1));
% scatter3(projectedC1_1(:,1),projectedC1_1(:,2),projectedC1_1(:,3),'r')% scatter3(projectedC2_1(:,1),projectedC2_1(:,2),projectedC2_1(:,3),'b')
vector(1,:)=normalize(projectedC1_1-pointCoord(1,:));vector(2,:)=normalize(projectedC2_1-pointCoord(1,:));
angle(1)=computeAngle(projectedC1_1,pointCoord(1,:),projectedC2_1);angle(2)=computeAngle(projectedC2_1,pointCoord(1,:),projectedC1_1);
%%%%%%%%%%%%%%generate circle point for branch 1
for i=1:2
    for h=1:nWalls/2;
        teta=(h-1)*(2*angle(i)/nWalls);
        circlePoints1{i}(h,1:3)=pointCoord(1,:)+(cos(teta)*radius(1)*vector(i,:))+(radius(1)*sin(teta)*normalize(cross(pointTangent(1,:),vector(i,:))));
    end
    %     scatter3(circlePoints1{i}(:,1),circlePoints1{i}(:,2),circlePoints1{i}(:,3),'*k')
end
bif.circlePoints{1}=[circlePoints1{1};circlePoints1{2}];

%%%%%%%%%%%%%%generate circle point for branch 2
projectedC1_2=computeStartPoint(controlPoints(1,:),bifCoord,bifTangent(2,:),pointCoord(2,:),pointTangent(2,:),radius(2));
projectedC2_2=computeStartPoint(controlPoints(2,:),bifCoord,bifTangent(2,:),pointCoord(2,:),pointTangent(2,:),radius(2));
vector(1,:)=normalize(projectedC1_2-pointCoord(2,:));
vector(2,:)=normalize(projectedC2_2-pointCoord(2,:));
angle(1)=computeAngle(projectedC1_2,pointCoord(2,:),projectedC2_2);
angle(2)=computeAngle(projectedC2_2,pointCoord(2,:),projectedC1_2);
for i=1:2
    for h=1:nWalls/2
        teta=-(h-1)*(2*angle(i)/nWalls);
        circlePoints2{i}(h,1:3)=pointCoord(2,:)+(cos(teta)*radius(2)*vector(i,:))+(radius(2)*sin(teta)*normalize(cross(pointTangent(2,:),vector(i,:))));
    end
    %     scatter3(circlePoints2{i}(:,1),circlePoints2{i}(:,2),circlePoints2{i}(:,3),'*b')
end
bif.circlePoints{2}=[circlePoints2{1};circlePoints2{2}];

%%%%%%%%%%%%%%generate circle point for branch 1
projectedC1_3=computeStartPoint(controlPoints(1,:),bifCoord,bifTangent(3,:),pointCoord(3,:),pointTangent(3,:),radius(3));
projectedC2_3=computeStartPoint(controlPoints(2,:),bifCoord,bifTangent(3,:),pointCoord(3,:),pointTangent(3,:),radius(3));
vector(1,:)=normalize(projectedC1_3-pointCoord(3,:));
vector(2,:)=normalize(projectedC2_3-pointCoord(3,:));
angle(1)=computeAngle(projectedC1_3,pointCoord(3,:),projectedC2_3);
angle(2)=computeAngle(projectedC2_3,pointCoord(3,:),projectedC1_3);
for i=1:2
    for h=1:nWalls/2
        teta=(h-1)*(2*angle(i)/nWalls);
        circlePoints3{i}(h,1:3)=pointCoord(3,:)+(cos(teta)*radius(3)*vector(i,:))+(radius(3)*sin(teta)*normalize(cross(pointTangent(3,:),vector(i,:))));
    end
    %     scatter3(circlePoints3{i}(:,1),circlePoints3{i}(:,2),circlePoints3{i}(:,3),'*r')
end
bif.circlePoints{3}=[circlePoints3{1};circlePoints3{2}];
Temporary_Circle=zeros(9,(nWalls*3+1));
trans=bif.circlePoints{1}';
cir1=trans(:);
Temporary_Circle(1,:)=[spline(1) cir1'];
trans=bif.circlePoints{2}';
cir2=trans(:);
Temporary_Circle(4,:)=[spline(2) cir2'];
trans=bif.circlePoints{3}';
cir3=trans(:);
Temporary_Circle(7,:)=[spline(3) cir3'];
bif.circle=[previousCir;Temporary_Circle];

function [ptCoordMx, faceMx, bif,  faceMxSize, ptCoordMxSize]= computePtCoordMx(bif,~,bifurSitPoints,~,nWalls,nSteps,c0,c1,pointCoord,bifCoord,ptCoordMx, faceMx,spline,bifFaceMx,splineTerm,ptCoordMxSize,faceMxSize,wallDensity,fluidDensity,ifMoreThan90)
tempBifMx =[];
wall = 0.8;
fluid = 0.45;
index = [];
d = (nWalls*(wallDensity+fluidDensity))+((nWalls/4)+1)^2;
for k = 1:3 
    counter = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%parallel bezier curves generation

    for i = 1:nWalls;
        %%% find all parallel control points
        counter = counter+1;
        %         scatter3(c0_off(:,1),c0_off(:,2),c0_off(:,3),'*m');%         scatter3(c1_off(:,1),c1_off(:,2),c1_off(:,3),'*g');
        
        %         P0new = minDistPointBezier2(pointCoord(k,:),c0(k,:),c1(k,:),bifCoord,bif.circlePoints{k}(i,:));
        P1new = minDistPointBezier2(pointCoord(k,:),c0(k,:),c1(k,:),bifCoord,bifurSitPoints{k}(i,:));
        distT = pdist([pointCoord(k,:);P1new]);
        t_P1new = distT/(distT+pdist([P1new;bifCoord]));
        %         scatter3(P0new(:,1),P0new(:,2),P0new(:,3),'*k')%         scatter3(P1new(:,1),P1new(:,2),P1new(:,3),'*b')
        %%%% made shortened bezier curve
        [curve1p0,curve1c0,curve1c1,curve1p1,splitPoint] = BezierSubdevision(pointCoord(k,:),c0(k,:),c1(k,:),bifCoord,t_P1new);
        p0_new = curve1p0(1 ,:);
        c0_new = curve1c0(1 ,:);
        c1_new = curve1c1(1 ,:);
        p1_new = curve1p1(1 ,:);
        %%%% compute parallel bezier curve
        offsetDirC0_1 = bif.circlePoints{k}(i,:)- pointCoord(k,:);
        offsetDirC1_1 = bifurSitPoints{k}(i,:) - p1_new;
        offsetDirC0 = ((2/3)*offsetDirC0_1+(1/3)*offsetDirC1_1);
        offsetDirC1 = ((2/3)*offsetDirC1_1+(1/3)*offsetDirC0_1);
        C0_Parallel = c0_new + pdist([bifurSitPoints{k}(i,:);p1_new])*normalize(offsetDirC0);
        C1_Parallel = c1_new + pdist([bif.circlePoints{k}(i,:);pointCoord(k,:)])*normalize(offsetDirC1);
        %         scatter3(c0(:,1),c0(:,2),c0(:,3),'filled','m');%         scatter3(c1(:,1),c1(:,2),c1(:,3),'filled','g');
        %                  scatter3(C0_Parallel(:,1),C0_Parallel(:,2),C0_Parallel(:,3),'filled','m');%                  scatter3(C1_Parallel(:,1),C1_Parallel(:,2),C1_Parallel(:,3),'filled','g');
        c0_off{k}(i,:) = C0_Parallel;
        c1_off{k}(i,:) = C1_Parallel;
    end
end

%%%%%%%%%%%%%%%%%%%trasfereing bifurcation site points too keep biurcation
%%%%%%%%%%%%%%%%%%%continuity between branches

for j = 1:3
    if j ==1 && ifMoreThan90(j) == 1;
        for f=2:(nWalls/2)
            %%%if branch1-2 is over 90 degreee
            
            %%%%giving weigth to the mean
            dist1 = pdist([c1_off{1}(f,:);bifurSitPoints{1}(f,:)]);
            dist2 = pdist([c1_off{2}(f,:);bifurSitPoints{1}(f,:)]);
            bifsiteT= (dist2/(dist2+dist1))*c1_off{1}(f,:) + (dist1/(dist2+dist1))*c1_off{2}(f,:);
            %               bifsiteT= c1_off{1}(f,:) + c1_off{2}(f,:);
            %                 scatter3(c1_off{1}(f,1),c1_off{1}(f,2),c1_off{1}(f,3),'*b');scatter3(c1_off{2}(f,1),c1_off{2}(f,2),c1_off{2}(f,3),'*g');scatter3(bifurSitPoints{1}(f,1),bifurSitPoints{1}(f,2),bifurSitPoints{1}(f,3),'*m');scatter3(bifsiteT(:,1),bifsiteT(:,2),bifsiteT(:,3),'filled','r');
            BifPointTrans{j}(f,:)= bifsiteT;
            if f ==(nWalls/2)
                bifurSitPoints{1}(2:(nWalls/2),:) = BifPointTrans{1}(2:(nWalls/2),:);
                bifurSitPoints{2}(2:(nWalls/2),:) = BifPointTrans{1}(2:(nWalls/2),:);
                
            end
        end
    end
    if j == 2 && ifMoreThan90(j) == 1
        for f = 2:nWalls/2;
            c1_temp3 = flipud(c1_off{3});
            c1_temp2= flipud(c1_off{2});
            %              scatter3(c1_temp3(f,1),c1_temp3(f,2),c1_temp3(f,3),'*b');%              scatter3(c1_temp2(f,1),c1_temp2(f,2),c1_temp2(f,3),'*m');
            dist1 = pdist([c1_temp2(f,:);bifurSitPoints{2}(f,:)]);
            dist2 = pdist([c1_temp3(f,:);bifurSitPoints{2}(f,:)]);
            bifsiteT= (dist2/(dist2+dist1))*c1_temp2(f,:)+(dist1/(dist2+dist1))*c1_temp3(f,:);
            %              bifsiteT= c1_temp2(f,:)+c1_temp3(f,:); %             scatter3(bifsiteT(:,1),bifsiteT(:,2),bifsiteT(:,3),'filled','r');
            BifPointTrans{j}(f,:)= bifsiteT;
            if f ==(nWalls/2)
                bifurSitPoints{2}((nWalls/2)+2:end,:) = flipud(BifPointTrans{2}(2:(nWalls/2),:));
                bifurSitPoints{3}((nWalls/2)+2:end,:) = flipud(BifPointTrans{2}(2:(nWalls/2),:));
            end
        end
    end
    %         %% if branch 1-3 is over 90 degree
    if j == 3 && ifMoreThan90(j) == 1
        for f=2:nWalls/2
            c1_temp = flipud(c1_off{1});
            dist1 = pdist([c1_temp(f,:);bifurSitPoints{3}(f,:)]);
            dist2 = pdist([c1_off{3}(f,:);bifurSitPoints{3}(f,:)]);
            bifsiteT= (dist2/(dist2+dist1))*c1_temp(f,:)+(dist1/(dist2+dist1))*c1_off{3}(f,:);
            BifPointTrans{j}(f,:)= bifsiteT;
            if f ==(nWalls/2)
                bifurSitPoints{1}((nWalls/2)+2:end,:) = flipud(BifPointTrans{3}(2:(nWalls/2),:));
                bifurSitPoints{3}(2:(nWalls/2),:) = BifPointTrans{3}(2:(nWalls/2),:);
            end
            
        end
    end
    
end
axis equal

% %%%%update control points
contro1UsedPoints2 = [bifurSitPoints{1}(2,:);bifurSitPoints{1}(end,:);bifurSitPoints{3}(end,:)];
% scatter3(contro1UsedPoints2(:,1),contro1UsedPoints2(:,2),contro1UsedPoints2(:,3),'filled','g');
contro1UsedPoints1 = [bifurSitPoints{1}((nWalls/2),:);bifurSitPoints{1}((nWalls/2)+2,:);bifurSitPoints{3}((nWalls/2)+2,:)];
% scatter3(contro1UsedPoints1(:,1),contro1UsedPoints1(:,2),contro1UsedPoints1(:,3),'filled','m');
contorlPoint2New = mean(contro1UsedPoints2);
contorlPoint1New = mean(contro1UsedPoints1);
% scatter3(contorlPoint2New(:,1),contorlPoint2New(:,2),contorlPoint2New(:,3),'filled','g');
% scatter3(contorlPoint1New(:,1),contorlPoint1New(:,2),contorlPoint1New(:,3),'filled','m');

bifurSitPoints{1}(1,:) = contorlPoint2New;
bifurSitPoints{1}((nWalls/2)+1,:) = contorlPoint1New;
bifurSitPoints{2}(1,:) = contorlPoint2New;
bifurSitPoints{2}((nWalls/2)+1,:) = contorlPoint1New;
bifurSitPoints{3}(1,:) = contorlPoint2New;
bifurSitPoints{3}((nWalls/2)+1,:) = contorlPoint1New;
%%%%%%%%%%%%%%%%%%%%%main cell generation for each cross section, you can
%%%%%%%%%%%%%%%%%%%%%visuzlize it here as well
for k = 1:3
    %%%%%%%
    counter = 0;
    for i = 1:nWalls
        counter = counter +1;
        wallBezierr = computeBezier(bif.circlePoints{k}(i,:),c0_off{k}(i,:),bifurSitPoints{k}(i,:),c1_off{k}(i,:),nSteps(k));
        %         plot3(wallBezierr(:,1),wallBezierr(:,2),wallBezierr(:,3),'linewidth',1,'color','b'); hold all;
        wallBezier{k}(:,1:3,counter) = wallBezierr;
        if counter>1
            %     fill(wallBezierr,wallBezierrLast,'b');
            pointPlotSurface = [wallBezierr;flipud(wallBezierrLast)];
            % patch3(pointPlotSurface(:,1),pointPlotSurface(:,2),pointPlotSurface(:,3))
            %             fill3(pointPlotSurface(:,1),pointPlotSurface(:,2),pointPlotSurface(:,3),rand([1 3]))
            %         'FaceAlpha', 0.9
        end
        wallBezierrLast = wallBezierr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ptCoordMx developmet%%%%%%%%%%%%%%%%%%%
        position = computeBezier(pointCoord(k,:), c0(k,:), bifCoord, c1(k,:),nSteps(k));
        
    end
    
    
    for j=1:nSteps(k)
        xx = reshape(wallBezier{k}(j,1, 1:nWalls),1,nWalls);
        yy = reshape(wallBezier{k}(j,2,1:nWalls),1,nWalls);
        zz = reshape(wallBezier{k}(j,3,1:nWalls),1,nWalls);
        %          plot3([xx xx(1)],[yy yy(1)],[zz zz(1)],'linewidth',0.5,'color','k');
        %%%%%%%%%%%%%%%inner mesh generation function
        [inMesh, squareP,ptCoordMx,indexMx,ptCoordMxSize] = computeInnermesh([xx' yy' zz'],position(j,:),wallDensity,fluidDensity,wall,fluid,nWalls,ptCoordMx,ptCoordMxSize);
        
        if j==1
            indexMx=indexMx';
            bif.index = [bif.index;[[spline(k);zeros(size(indexMx,1)-1,1)],indexMx]];
            % bif.index = [spline(k);
        end
        %     %%%%%%end of inner mesh generation function
        %%%%%%%%%%%%%%writing ptCoordMx
        nStepsIndex = j;
        branchIndex = k;
        csDensity = fluidDensity + wallDensity;
        [r,~] = find(bif.circle(:,1) == spline(k));
        if j ==1
            rightHandTest = ifRightHand(inMesh{1}(1,2:4),inMesh{1}(2,2:4),position(1,:),(c0(k,:)-position(1,:)));
        end
        if j ~= 1
            [faceMx, bif,bifFaceMx,faceMxSize] = inletFaceMaker (inMesh,faceMx, nStepsIndex, bif,nWalls,r,branchIndex,nSteps,tempBifMx,bifFaceMx,csDensity,squareP,d,ptCoordMx,rightHandTest,j,ptCoordMxSize,faceMxSize);
        end
        if j ~=nSteps(k)
            [faceMx, bif,bifFaceMx,faceMxSize] = innerFaceMxMaker (inMesh,faceMx, nStepsIndex, bif,nWalls,r,branchIndex,nSteps,tempBifMx,bifFaceMx,csDensity,squareP,d,ptCoordMx,rightHandTest,ptCoordMxSize,faceMxSize);
        end
        if j == nSteps(k)
            index  = computeBifIndex(index,inMesh,squareP,nWalls,csDensity,k);
            xxx{k} = reshape(wallBezier{k}(j+1,1, 1:nWalls),1,nWalls);
            yyy{k} = reshape(wallBezier{k}(j+1,2,1:nWalls),1,nWalls);
            zzz{k} = reshape(wallBezier{k}(j+1,3,1:nWalls),1,nWalls);
        end
    end
end

[inMesh, squareP,ptCoordMx,~,ptCoordMxSize] = computeInnermesh([xxx{1}' yyy{1}' zzz{1}'],position(j+1,:),wallDensity,fluidDensity,wall,fluid,nWalls,ptCoordMx,ptCoordMxSize);

indexTemp  = computeBifIndex(index,inMesh,squareP,nWalls,csDensity,1);
upP = index.L;middleP = indexTemp.L;downP = index.LR; branch = 1;
[faceMx, faceMxSize] = makeBifFace(upP,middleP,downP,faceMx,csDensity,nWalls,branch,ptCoordMxSize,faceMxSize,d);
upP = index.R;
middleP = indexTemp.R;
downP = index.RL;
branch = 2;
[faceMx, faceMxSize] = makeBifFace(upP,middleP,downP,faceMx,csDensity,nWalls,branch,ptCoordMxSize,faceMxSize,d);
[middleIndex, ptCoordMx,ptCoordMxSize] = computeInnermesh1([xxx{2}' yyy{2}' zzz{2}'],position(j+1,:),wallDensity,fluidDensity,wall,fluid,nWalls,ptCoordMx,inMesh, squareP,ptCoordMxSize);
downP = index.LL;
middleP = middleIndex;
upP = index.RR;
branch=3;
[faceMx, faceMxSize] = makeBifFace(upP,middleP,downP,faceMx,csDensity,nWalls,branch,ptCoordMxSize,faceMxSize,d);

function [middleIndex, ptCoordMx,ptCoordMxSize]=computeInnermesh1(points,center,wallDensity,fluidDensity,wall,fluid,nWalls,ptCoordMx,inMesh, squareP,ptCoordMxSize)
index =[];
direction = zeros(nWalls,3);
totalP = (nWalls*(wallDensity + fluidDensity))+((nWalls/4)+1)^2;
tempPtCoord =[];
for h=1:wallDensity+fluidDensity
    u2 = inMesh{1,h}(7*(nWalls/8)+1,:);
    index = [index;u2];
end
index = [index; squareP{1,(nWalls/8)+1}(1:end,:)];
criticalDir = zeros(8,3);
criticalP1 = zeros(8,3);
criticalP2 = zeros(8,3);
for h = wallDensity+fluidDensity:-1:1
    u1=inMesh{1,h}(3*(nWalls/8)+1,:);
    index=[index; u1 ];
end
d = 0:7;
criticalP = points((1+(nWalls/8)*d),:);%%%It should be 8

for y = 1:nWalls
    direction(y,:) = points(y,:) - center;
end
for y = 1:8
    criticalDir(y,:) = (criticalP(y,:) - center);
    criticalP1(y,:) = center + criticalDir(y,:).*wall;
    if mod(y,2) == 1
        criticalP2(y,:) = center + criticalDir(y,:).*fluid;
    else
        criticalP2(y,:) = center + criticalDir(y,:).*fluid*(4.2/3.5);
    end
end

%%%lumen DOMAIN MESH
for b = wallDensity+1:-1:1
    for u = 1:nWalls
        wallP{b}(u,:) = center + direction(u,:).*(wall+(((1-wall)/wallDensity)*(b-1)));
        %           scatter3(wallP{b}(u,1),wallP{b}(u,2),wallP{b}(u,3));
    end
end
%%%%line Point generation
counter = 0;
lineP = cell(1,nWalls/4);
for h = 1:2:8
    counter = counter+1;
    DirLine = criticalP2(h,:)-center;
    for u = 1:(nWalls/8)
        lineP{counter}(u,:) = center + 8*DirLine*u/nWalls;
        %             scatter3(lineP{counter}(u,1), lineP{counter}(u,2),lineP{counter}(u,3));
    end
end

cp2 =  [criticalP2;criticalP2(1,:)];
dir = zeros(8,3);
circle = cell(1,8);
square = cell(1,8);
dirSC = cell(1,8);

for b = 1:8
    dir(b,:) = cp2(b+1,:) - cp2(b,:);
    circle{b}(1:(nWalls/8),:) =  wallP{1}((1+(b-1)*(nWalls/8)):b*(nWalls/8),:);
    for u = 1:(nWalls/8)
        
        square{b}(u,:) = cp2(b,:) + dir(b,:)*(8*(u-1)/nWalls);     %%dir(b,:).*(dist(b)*(u-1));
        %           scatter3(square{b}(u,1),square{b}(u,2),square{b}(u,3))
        
        dirSC{b}(u,:) = circle{b}(u,:) - square{b}(u,:);
        
    end
end


fluidP = cell(1,8);
for h=1:8
    for u = 1:nWalls/8
        for c = 1:(fluidDensity+1)
            fluidP{h}(u,:,c) = square{h}(u,:) + dirSC{h}(u,:)*((c-1)/fluidDensity);
            %                  scatter3(fluidP{h}(u,1,c),fluidP{h}(u,2,c),fluidP{h}(u,3,c))
        end
    end
end

%%%%%Indexing for the square inside the circle 
%square 1
terminal = [flipud(lineP{2});center];
start = [square{2}(1,:);flipud(square{1})];
dir = terminal-start;
sPart1 = cell(1,(nWalls/8)+1);
sPart2 = cell(1,(nWalls/8)+1);
sPart4 = cell(1,(nWalls/8)+1);
sPart3 = cell(1,(nWalls/8)+1);
for k = 1:(nWalls/8)+1
    for u =1:(nWalls/8)+1
        sPart1{k}(u,:) = start(k,:)+ dir(k,:)*(u-1)/(nWalls/8);
    end
end
%%mini square2
terminal=[square{4};square{5}(1,:)];
start = [flipud(lineP{2});center];
dir = terminal-start;
for k = 1:(nWalls/8)+1
    for u = 1:(nWalls/8)
        sPart2{k}(u,:) = start(k,:) + dir(k,:)*(u)/(nWalls/8);
        
    end
end
%%mini square4
start = flipud(square{8});
terminal = lineP{4};
dir = terminal-start;
for k = 1:(nWalls/8)
    for u = 1:(nWalls/8)+1
        sPart4{k}(u,:) = start(k,:) + dir(k,:)*(u-1)/(nWalls/8);
        
    end
end

%%mini square3
start = lineP{4};
terminal = [square{5}(2:end,:);square{6}(1,:)];
dir = terminal-start;
for k = 1:(nWalls/8)
    for u =1:(nWalls/8)
        sPart3{k}(u,:) = start(k,:)+ dir(k,:)*(u)/(nWalls/8);
        
    end
end

for k=1:(nWalls/8)+1
    squareP{k} = [sPart1{k};sPart2{k}];
    % scatter3(squareP{k}(:,1),squareP{k}(:,2),squareP{k}(:,3))
end
for k= 1:(nWalls/8)
    squareP{k+((nWalls/8)+1)} = [sPart4{k};sPart3{k}];
    % scatter3(squareP{k+((nWalls/8)+1)}(:,1),squareP{k+((nWalls/8)+1)}(:,2),squareP{k+((nWalls/8)+1)}(:,3))
end





%%%csWritting

for k = 1:wallDensity
    inMesh{k+1} = wallP{(wallDensity+1)-k};
end
for k = (fluidDensity+1):-1:1
    f=[];
    for u = 1:8
        b = fluidP{u}(1:(nWalls/8),1:3,k);
        f = [f;b];
    end
    inMesh{fluidDensity+2+wallDensity-k} = f;

end
inMesh{1} = points;
tempSize = ptCoordMxSize;
for f = 1:(wallDensity+fluidDensity)
    ptIndex = ((tempSize+1):(tempSize+((nWalls/2)-1)))';
    b=flipud(inMesh{f}(((nWalls/2)+2):nWalls, 1:3));
    tempPtCoord = [tempPtCoord;b];
    middleIndex(:,f) = [index(f,1);ptIndex;index(size(index,1)-f+1,1)]; %#ok<AGROW>
    tempSize = tempSize +(nWalls/2)-1;
end




counter=0;
for f = (nWalls/4+1):-1:(nWalls/8+2)
    counter = counter +1;
    ptIndex = ((tempSize+1):(tempSize+(nWalls/4)+1))';
    tempPtCoord = [tempPtCoord; squareP{f}];
    tempSize = tempSize+(nWalls/4)+1;
    squareP{f} = [ptIndex,squareP{f}(:,1:3)];
    middleIndex(1:length(ptIndex),counter+wallDensity+fluidDensity) = ptIndex;
end
ptCoordMx(ptCoordMxSize+1:(ptCoordMxSize+length(tempPtCoord)),:)= tempPtCoord;
ptCoordMxSize = ptCoordMxSize+length(tempPtCoord);
middleIndex(1:length(ptIndex),end+1) = flipud(index(wallDensity+fluidDensity+1,1):index(end-(wallDensity+fluidDensity),1));
middleIndex(:,wallDensity+fluidDensity+1:end)=  fliplr(middleIndex(:,wallDensity+fluidDensity+1:end));
function [faceMx, faceMxSize] = makeBifFace(upP,middleP,downP,faceMx,csDensity,nWalls,branch,ptCoordMxSize,faceMxSize,d)
iFace1 = faceMxSize.in+1;
nS1 = faceMxSize.sur+1;
%%%uP,mP,dP are points and cells around the square are, the first row is
%%%the points and the second raw are the cell indices.
w = nWalls/4;
uP(1,:) = [upP(1,csDensity+1:end),upP(2:w,end)',fliplr(upP(w+1,csDensity+1:end))];
if branch ==1
    uP(2,:) = [upP(1,csDensity+2:end),upP(1:w,end)',fliplr(upP(w,csDensity+2:end)),0];
elseif branch == 2 || branch == 3
    uP(2,:) = [upP(1,csDensity+1:end-1),upP(1:w,end-1)',fliplr(upP(w,csDensity+1:end-1)),0];
end
mP = [middleP(1,csDensity+1:end),middleP(2:w,end)',fliplr(middleP(w+1,csDensity+1:end))];

dP(1,:) = [downP(1,csDensity+1:end),downP(2:w,end)',fliplr(downP(w+1,csDensity+1:end))];
if branch ==1 || branch == 2
    dP(2,:) = [downP(1,csDensity+2:end),downP(1:w,end)',fliplr(downP(w,csDensity+2:end)),0];
else
    dP(2,:) =  [downP(1,csDensity+1:end-1),downP(1:w,end-1)',fliplr(downP(w,csDensity+1:end-1)),0];
end
localFaceMxIn = zeros(d*2,6);
localFaceMxSur = zeros(nWalls,6);
nS = 1;
iFace = 1;
if branch ==1
    for t = 1:csDensity
        i = 0;
        for a = 1:(nWalls/2)
            if t == 1
                i = i+1;
                localFaceMxIn(iFace,:) =[middleP(i,t) middleP(i,t+1) middleP(i+1,t+1) middleP(i+1,t) upP(i,t) downP(i,t)];
                localFaceMxSur(nS,:) =[middleP(i,t) middleP(i+1,t) upP(i+1,t) upP(i,t) upP(i,t) 0];
                localFaceMxSur(nS+1,:) =[middleP(i,t) downP(i,t) downP(i+1,t) middleP(i+1,t)   downP(i,t) 0 ];
                nS = nS+2;
                iFace = iFace+1;
                if a ==1
                    localFaceMxIn(iFace,:) = [middleP(i,t) middleP(i,t+1) downP(i,t+1) downP(i,t)  downP(i,t) downP(i,t)-1];
                    iFace = iFace+1;
                else
                    localFaceMxIn(iFace,:) = [middleP(i,t) middleP(i,t+1) upP(i,t+1) upP(i,t) upP(i-1,t) upP(i,t)];
                    localFaceMxIn(iFace+1,:)=[middleP(i,t) middleP(i,t+1) downP(i,t+1) downP(i,t)  downP(i,t) downP(i-1,t)];
                    iFace = iFace+2;
                end
                if a ==(nWalls/2)
                    localFaceMxIn(iFace,:) = [middleP(i+1,t) middleP(i+1,t+1) downP(i+1,t+1) downP(i+1,t)  downP(i+1,t) downP(i,t)];
                    iFace = iFace+1;
                end
            elseif (t > 1 && t < csDensity)
                i = i+1;%%%%% irad darad barrasi shavad
                localFaceMxIn(iFace,:) = [middleP(i,t) middleP(i,t+1)  middleP(i+1,t+1) middleP(i+1,t) upP(i,t) downP(i,t)];
                localFaceMxIn(iFace+1,:) = [middleP(i,t) middleP(i+1,t) upP(i+1,t) upP(i,t) upP(i,t) upP(i,t-1) ];
                localFaceMxIn(iFace+2,:) = [middleP(i,t) downP(i,t) downP(i+1,t) middleP(i+1,t) downP(i,t) downP(i,t-1) ];
                iFace = iFace+3;
                if a ==1
                    localFaceMxIn(iFace,:) = [middleP(i,t) middleP(i,t+1) downP(i,t+1) downP(i,t)  downP(i,t) downP(i,t)-1];
                    iFace = iFace+1;
                else
                    localFaceMxIn(iFace,:) = [middleP(i,t) middleP(i,t+1) upP(i,t+1) upP(i,t) upP(i-1,t) upP(i,t)];
                    localFaceMxIn(iFace+1,:) = [middleP(i,t) middleP(i,t+1) downP(i,t+1) downP(i,t)  downP(i,t) downP(i-1,t)];
                    iFace = iFace+2;
                end
                if a ==(nWalls/2)
                    localFaceMxIn(iFace,:) = [middleP(i+1,t) middleP(i+1,t+1) downP(i+1,t+1) downP(i+1,t)  downP(i+1,t) downP(i,t)];
                    iFace = iFace+1;
                end
                
                %%%cells between circle and square
                
            elseif t == csDensity
                i = i+1;
                localFaceMxIn(iFace,:) = [middleP(i,t) mP(1,i) mP(1,i+1) middleP(i+1,t) upP(i,t) downP(i,t)];
                localFaceMxIn(iFace+1,:) = [middleP(i,t) middleP(i+1,t) upP(i+1,t) upP(i,t) upP(i,t) upP(i,t-1)];
                localFaceMxIn(iFace+2,:) = [middleP(i,t) downP(i,t) downP(i+1,t) middleP(i+1,t) downP(i,t) downP(i,t-1)];
                localFaceMxIn(iFace+3,:) = [mP(1,i) uP(1,i) uP(1,i+1) mP(1,i+1) upP(i,t) uP(2,i)];
                localFaceMxIn(iFace+4,:) = [mP(1,i) mP(1,i+1) dP(1,i+1) dP(1,i) downP(i,t)  dP(2,i)];
                iFace = iFace+5;
                if a == 1
                    localFaceMxIn(iFace,:) = [middleP(i,t) mP(1,i) dP(1,i) downP(i,t)  downP(i,t) downP(i,t)-1];
                    iFace = iFace+1;
                else
                    localFaceMxIn(iFace,:) = [middleP(i,t) upP(i,t) uP(1,i) mP(1,i) upP(i,t) upP(i-1,t)];
                    localFaceMxIn(iFace+1,:) = [middleP(i,t) mP(1,i) dP(1,i) downP(i,t)  downP(i,t) downP(i-1,t)];
                    iFace = iFace+2;
                end
                if a ==(nWalls/2)
                    localFaceMxIn(iFace,:) = [middleP(i+1,t) mP(1,i+1) dP(1,i+1) downP(i+1,t)  downP(i+1,t) downP(i,t)];
                    iFace = iFace+1;
                end
            end
        end
    end
    
    
    for t = (csDensity+1):(csDensity+(nWalls/8))
        i = 0;
        for a = 1:(nWalls/4);
            if a == (nWalls/4)
                i = i+1;
                localFaceMxIn(iFace,:) = [middleP(i,t) middleP(i+1,t) middleP(i+1,t+1) middleP(i,t+1) upP(i,t+1) downP(i,t+1)];
                localFaceMxIn(iFace+1,:) = [middleP(i,t) middleP(i+1,t) downP(i+1,t)  downP(i,t) downP(i,t+1) downP(i,t)];
                iFace = iFace+2;
                if t ~= csDensity+1
                    localFaceMxIn(iFace,:) = [middleP(i,t) upP(i,t) upP(i+1,t) middleP(i+1,t)  upP(i,t+1) upP(i,t)]; %%%% barrasi shavd
                    %        faceMx.in(iFace+1,:) = [middleP(i,t) middleP(i+1,t) downP(i+1,t)  downP(i,t) downP(i,t+1) downP(i,t)];%%%% barrasi shavd
                    iFace = iFace+1;
                end
            else
                i = i+1;
                localFaceMxIn(iFace,:) = [middleP(i,t) middleP(i+1,t) middleP(i+1,t+1) middleP(i,t+1) upP(i,t+1) downP(i,t+1)];
                localFaceMxIn(iFace+1,:) = [middleP(i+1,t) upP(i+1,t) upP(i+1,t+1) middleP(i+1,t+1) upP(i,t+1) upP(i+1,t+1)];
                localFaceMxIn(iFace+2,:) = [middleP(i+1,t) middleP(i+1,t+1) downP(i+1,t+1) downP(i+1,t) downP(i,t+1) downP(i+1,t+1)];
                localFaceMxIn(iFace+3,:) = [middleP(i,t) middleP(i+1,t) downP(i+1,t)  downP(i,t) downP(i,t+1) downP(i,t)];
                iFace = iFace+4;
                if t ~= csDensity+1
                    localFaceMxIn(iFace,:) = [middleP(i,t) upP(i,t) upP(i+1,t) middleP(i+1,t)  upP(i,t+1) upP(i,t)]; %%%% barrasi shavd
                    iFace = iFace+1;
                end
            end
            
        end
    end
end


if branch == 2
    for t = 1:csDensity
        i = 0;
        for a = 1:(nWalls/2)
            if t == 1
                i = i+1;
                localFaceMxIn(iFace,:) = [middleP(i,t) middleP(i,t+1) middleP(i+1,t+1) middleP(i+1,t) downP(i,t) upP(i+1,t)];
                localFaceMxIn(iFace+1,:) = [middleP(i,t) middleP(i,t+1) upP(i,t+1) upP(i,t) upP(i+1,t) upP(i,t)];
                iFace = iFace+2;
                if a ~= 1
                    localFaceMxIn(iFace,:) = [middleP(i,t) middleP(i,t+1) downP(i,t+1) downP(i,t) downP(i-1,t) downP(i,t) ];
                    iFace = iFace+1;
                end
                localFaceMxSur(nS,:) = [middleP(i,t) middleP(i+1,t) upP(i+1,t) upP(i,t) 0 upP(i+1,t)];
                localFaceMxSur(nS+1,:) = [middleP(i,t) downP(i,t) downP(i+1,t) middleP(i+1,t) 0 downP(i,t)];
                
                nS = nS+2;
                if a == (nWalls/2)
                    localFaceMxIn(iFace,:) = [middleP(i+1,t) middleP(i+1,t+1) upP(i+1,t+1) upP(i+1,t) upP(i+1,t)-1 upP(i+1,t)];
                    iFace = iFace+1;
                end
            elseif (t > 1 && t < csDensity)
                %%%IRAD DARE
                i = i+1;
                localFaceMxIn(iFace,:) = [middleP(i,t) middleP(i,t+1) middleP(i+1,t+1) middleP(i+1,t) downP(i,t) upP(i+1,t)];
                localFaceMxIn(iFace+1,:) = [middleP(i,t) middleP(i,t+1) upP(i,t+1) upP(i,t) upP(i+1,t) upP(i,t)];%
                localFaceMxIn(iFace+2,:) = [middleP(i,t) middleP(i+1,t) upP(i+1,t) upP(i,t)  upP(i+1,t-1) upP(i+1,t)]; %%barrasi shavad
                localFaceMxIn(iFace+3,:) = [middleP(i,t) downP(i,t) downP(i+1,t) middleP(i+1,t)  downP(i,t-1) downP(i,t)]; %%barrasi shavad
                if a == 1
                    iFace = iFace+4;
                else
                    localFaceMxIn(iFace+4,:) = [middleP(i,t) middleP(i,t+1) downP(i,t+1) downP(i,t) downP(i-1,t) downP(i,t) ];
                    iFace = iFace+5;
                end
                if a == (nWalls/2)
                    localFaceMxIn(iFace,:) = [middleP(i+1,t) middleP(i+1,t+1) upP(i+1,t+1) upP(i+1,t) upP(i+1,t)-1 upP(i+1,t)];
                    iFace = iFace+1;
                end
                
                
                %%%cells between circle and square
                
            elseif t == csDensity
                i = i+1;
                localFaceMxIn(iFace,:) = [middleP(i,t) mP(1,i) mP(1,i+1) middleP(i+1,t)  downP(i,t) upP(i+1,t)];
                localFaceMxIn(iFace+1,:) = [middleP(i,t) upP(i,t) uP(1,i) mP(1,i) upP(i,t) upP(i+1,t)];
                localFaceMxIn(iFace+2,:) = [middleP(i,t) middleP(i+1,t) upP(i+1,t) upP(i,t)  upP(i+1,t-1) upP(i+1,t)];%%irad
                localFaceMxIn(iFace+3,:) = [middleP(i,t) downP(i,t) downP(i+1,t) middleP(i+1,t)  downP(i,t-1) downP(i,t)];
                localFaceMxIn(iFace+4,:) = [mP(1,i) uP(1,i) uP(1,i+1) mP(1,i+1) uP(2,i) upP(i+1,t)];
                localFaceMxIn(iFace+5,:) = [mP(1,i) mP(1,i+1) dP(1,i+1) dP(1,i) dP(2,i) downP(i,t)];
                if a ==1
                    iFace = iFace+6;
                else
                    localFaceMxIn(iFace+6,:) = [middleP(i,t) mP(1,i) dP(1,i) downP(i,t)  downP(i-1,t) downP(i,t)];
                    iFace = iFace+7;
                end
                
                if a == (nWalls/2)
                    localFaceMxIn(iFace,:)= [middleP(i+1,t)  mP(1,i+1) uP(1,i+1) upP(i+1,t) upP(i+1,t)-1 upP(i+1,t)];
                    iFace = iFace+1;
                end
                
                
                
            end
        end
    end
    
    for t = (csDensity+1):(csDensity+(nWalls/8))
        i = 0;
        for a = 1:(nWalls/4);
            if a == (nWalls/4)
                i = i+1;
                localFaceMxIn(iFace,:) = [middleP(i,t) middleP(i,t)+1 middleP(i,t+1)+1 middleP(i,t+1) downP(i,t+1)  upP(i,t)];
                localFaceMxIn(iFace+1,:) = [middleP(i,t) upP(i,t) upP(i,t)+1 middleP(i,t)+1  upP(i,t)-((nWalls/4)+1) upP(i,t)];
                iFace = iFace+2;
                if t ~= csDensity+1
                    localFaceMxIn(iFace,:) = [middleP(i,t) middleP(i,t)+1 downP(i,t)+1  downP(i,t) downP(i,t) downP(i,t+1)];
                    iFace = iFace+1;
                end
            else
                i = i+1;
                localFaceMxIn(iFace,:) = [middleP(i,t) middleP(i,t)+1 middleP(i,t+1)+1 middleP(i,t+1) downP(i,t+1) upP(i,t) ];
                localFaceMxIn(iFace+1,:) = [middleP(i,t) upP(i,t) upP(i,t)+1 middleP(i,t)+1   upP(i,t)-((nWalls/4)+1) upP(i,t)];
                localFaceMxIn(iFace+2,:) = [middleP(i,t)+1 upP(i,t)+1 upP(i,t+1)+1 middleP(i,t+1)+1 upP(i+1,t) upP(i,t)];
                localFaceMxIn(iFace+3,:) = [middleP(i,t)+1 middleP(i,t+1)+1 downP(i,t+1)+1 downP(i,t)+1  downP(i+1,t+1) downP(i,t+1)];
                iFace = iFace+4;
                if t ~= csDensity+1
                    localFaceMxIn(iFace,:) = [middleP(i,t) middleP(i,t)+1 downP(i,t)+1  downP(i,t)  downP(i,t) downP(i,t+1)];
                    iFace = iFace+1;
                end
            end
            
        end
    end
end
if branch == 3
    for t=1:csDensity
        i=0;
        for a = 1:(nWalls/2)
            if t == 1
                i=i+1;
                localFaceMxIn(iFace,:) = [middleP(i,t) middleP(i+1,t) middleP(i+1,t+1) middleP(i,t+1) upP(i+1,t) downP(i+1,t)];%%barrasi shavad
                localFaceMxIn(iFace+1,:) = [middleP(i,t) middleP(i,t+1) upP(i,t+1) upP(i,t) upP(i+1,t) upP(i,t)];
                localFaceMxSur(nS,:) = [middleP(i,t) middleP(i+1,t) upP(i+1,t) upP(i,t)  0 upP(i+1,t)];
                localFaceMxSur(nS+1,:) = [middleP(i,t) downP(i,t) downP(i+1,t) middleP(i+1,t) 0 downP(i+1,t) ];
                nS = nS+2;
                iFace = iFace +2;
                if a ~= 1
                    localFaceMxIn(iFace,:) = [middleP(i,t) middleP(i,t+1) downP(i,t+1) downP(i,t) downP(i,t) downP(i+1,t)];
                    iFace = iFace +1;
                end
                if a == (nWalls/2)
                    localFaceMxIn(iFace,:) =[middleP(i+1,t) middleP(i+1,t+1) upP(i+1,t+1) upP(i+1,t)  upP(i+1,t)-1 upP(i+1,t)];
                    iFace = iFace +1;
                end
                
                
                
            elseif (t > 1 && t < csDensity)
                i=i+1;
                localFaceMxIn(iFace,:) = [middleP(i,t) middleP(i+1,t) middleP(i+1,t+1) middleP(i,t+1) upP(i+1,t) downP(i+1,t)];
                localFaceMxIn(iFace+1,:) = [middleP(i,t) middleP(i,t+1) upP(i,t+1) upP(i,t) upP(i+1,t) upP(i,t)];
                localFaceMxIn(iFace+2,:) = [middleP(i,t) middleP(i+1,t) upP(i+1,t) upP(i,t)  upP(i+1,t-1) upP(i+1,t)];
                localFaceMxIn(iFace+3,:) = [middleP(i,t) downP(i,t) downP(i+1,t) middleP(i+1,t)  downP(i+1,t-1) downP(i+1,t)];
                iFace = iFace +4;
                
                if a ~= 1
                    localFaceMxIn(iFace,:) = [middleP(i,t) middleP(i,t+1) downP(i,t+1) downP(i,t) downP(i,t) downP(i+1,t)];
                    iFace = iFace +1;
                end
                if a == (nWalls/2)
                    localFaceMxIn(iFace,:) =[middleP(i+1,t) middleP(i+1,t+1) upP(i+1,t+1) upP(i+1,t)    upP(i+1,t)-1 upP(i+1,t)];
                    iFace = iFace +1;
                end
            elseif  t == csDensity
                i = i+1;
                localFaceMxIn(iFace,:) = [middleP(i,t) middleP(i+1,t) mP(i+1) mP(i) upP(i+1,t) downP(i+1,t)];
                localFaceMxIn(iFace+1,:) = [middleP(i,t) mP(i) dP(1,i) downP(i,t) downP(i,t) downP(i+1,t) ];
                localFaceMxIn(iFace+2,:) = [mP(i) mP(i+1) uP(1,i+1) uP(1,i) upP(i+1,t) uP(2,i) ];
                localFaceMxIn(iFace+3,:) = [mP(i) mP(i+1) dP(1,i+1) dP(1,i) dP(2,i) downP(i+1,t)];
                localFaceMxIn(iFace+4,:) = [middleP(i,t) middleP(i+1,t) upP(i+1,t) upP(i,t)  upP(i+1,t-1) upP(i+1,t)];
                localFaceMxIn(iFace+5,:) = [middleP(i,t) downP(i,t) downP(i+1,t) middleP(i+1,t)  downP(i+1,t-1) downP(i+1,t)];
                iFace = iFace +6;
                
                if a~= 1
                    localFaceMxIn(iFace,:) = [middleP(i,t) mP(i) uP(1,i) upP(i,t) upP(i+1,t) upP(i,t)];
                    iFace = iFace +1;
                end
                if a == (nWalls/2)
                    localFaceMxIn(iFace,:) =[middleP(i+1,t) mP(i+1) uP(1,i+1) upP(i+1,t)  upP(i+1,t)-1 upP(i+1,t)];
                    iFace = iFace +1;
                end
            end
        end
    end
    for t = (csDensity+1):(csDensity+(nWalls/8))
        i = 0;
        for a = 1:(nWalls/4);
            if a == (nWalls/4)
                i = i+1;
                localFaceMxIn(iFace,:) = [middleP(i,t) middleP(i+1,t) middleP(i+1,t+1) middleP(i,t+1)  downP(i,t) upP(i,t)];
                localFaceMxIn(iFace+1,:) = [middleP(i,t) upP(i,t) upP(i+1,t) middleP(i+1,t)  upP(i,t)-((nWalls/4)+1) upP(i,t)];
                iFace = iFace+2;
                if t ~= csDensity+1
                    %        faceMx.in(iFace,:) = [middleP(i,t) upP(i,t) upP(i+1,t) middleP(i+1,t)  upP(i,t-1) upP(i,t)]; %%%% barrasi shavd
                    localFaceMxIn(iFace,:) = [middleP(i,t) middleP(i+1,t) downP(i+1,t)  downP(i,t) downP(i,t-1) downP(i,t)];%%%% barrasi shavd
                    iFace = iFace+1;
                end
            else
                i = i+1;
                localFaceMxIn(iFace,:) = [middleP(i,t) middleP(i+1,t) middleP(i+1,t+1) middleP(i,t+1)  downP(i,t) upP(i,t)];
                localFaceMxIn(iFace+1,:) = [middleP(i,t) upP(i,t) upP(i+1,t) middleP(i+1,t)  upP(i,t)-((nWalls/4)+1) upP(i,t)];%%%i changed it!
                localFaceMxIn(iFace+2,:) = [middleP(i+1,t) upP(i+1,t) upP(i+1,t+1) middleP(i+1,t+1)  upP(i+1,t) upP(i,t)];
                localFaceMxIn(iFace+3,:) = [middleP(i+1,t) middleP(i+1,t+1) downP(i+1,t+1) downP(i+1,t)   downP(i+1,t) downP(i,t)];
                iFace = iFace+4;
                if t ~= csDensity+1
                    %        faceMx.in(iFace,:) = [middleP(i,t) upP(i,t) upP(i+1,t) middleP(i+1,t)  upP(i,t-1) upP(i,t)]; %%%% barrasi shavd
                    localFaceMxIn(iFace,:) = [middleP(i,t) middleP(i+1,t) downP(i+1,t)  downP(i,t) downP(i,t-1) downP(i,t)];%%%% barrasi shavd
                    iFace = iFace+1;
                end
            end
        end
    end
end
localFaceMxSur = localFaceMxSur(any(localFaceMxSur,2),:);
localFaceMxIn = localFaceMxIn(any(localFaceMxIn,2),:);
faceMx.in(faceMxSize.in+1:(faceMxSize.in+length(localFaceMxIn)),:) = localFaceMxIn;
faceMx.sur(faceMxSize.sur+1:(faceMxSize.sur+length(localFaceMxSur)),:) = localFaceMxSur;
faceMxSize.in = faceMxSize.in + length(localFaceMxIn);
faceMxSize.sur = faceMxSize.sur + length(localFaceMxSur);
function index  = computeBifIndex(index,inMesh,squareP,nWalls,csDensity,k)


if k == 1
    index.L = zeros((nWalls/2)+1,csDensity+(nWalls/8)+1);
    tempIndex = zeros(nWalls,4);
    for f = 1:(csDensity)
        tempIndex(:,f) = inMesh{f}(:,1);
    end
    tempIndex = [tempIndex(7*(nWalls/8)+1:end,:);tempIndex(1:(7*(nWalls/8))+1,:)];
    index.L = (tempIndex(1:(nWalls/2)+1,:));
    index.R = flipud((tempIndex((nWalls/2)+1:end,:)));
    
    for f= 1:(nWalls/4)+1
        tempIndex1(:,f) = squareP{f}(:,1);
    end
    for f= 1:(nWalls/8)+1
        index.L(1:(nWalls/4)+1,f+csDensity) = tempIndex1(:,(nWalls/8)+2-f);
        index.R(1:(nWalls/4)+1,f+csDensity) = tempIndex1(:,f+(nWalls/8));
    end
end


if k == 2
    for f = 1:(csDensity)
        tempIndex2(:,f) = inMesh{f}(:,1);
    end
    tempIndex2 = [tempIndex2(7*(nWalls/8)+1:end,:);tempIndex2(1:(7*(nWalls/8))+1,:)];
    index.LR = (tempIndex2(1:(nWalls/2)+1,:));
    index.LL = flipud((tempIndex2((nWalls/2)+1:end,:)));
    for f= 1:(nWalls/4)+1
        tempIndex3(:,f) = squareP{f}(:,1);
    end
    for f= 1:(nWalls/8)+1
        index.LR(1:(nWalls/4)+1,f+csDensity) = tempIndex3(:,(nWalls/8)+2-f);
        index.LL(1:(nWalls/4)+1,f+csDensity) = tempIndex3(:,f+(nWalls/8));
    end
end


if k == 3
    for f = 1:(csDensity)
        tempIndex4(:,f) = inMesh{f}(:,1);
    end
    tempIndex4 = [tempIndex4(7*(nWalls/8)+1:end,:);tempIndex4(1:(7*(nWalls/8))+1,:)];
    index.RL = (tempIndex4(1:(nWalls/2)+1,:));
    index.RR = flipud((tempIndex4((nWalls/2)+1:end,:)));
    for f= 1:(nWalls/4)+1
        tempIndex5(:,f) = squareP{f}(:,1);
    end
    for f= 1:(nWalls/8)+1
        index.RL(1:(nWalls/4)+1,f+csDensity) = tempIndex5(:,(nWalls/8)+2-f);
        index.RR(1:(nWalls/4)+1,f+csDensity) = tempIndex5(:,f+(nWalls/8));
    end
end

function [inMesh,squareP,ptCoordMx,indexMx,ptCoordMxSize] = computeInnermesh(points,center,wallDensity,fluidDensity,wall,fluid,nWalls,ptCoordMx,ptCoordMxSize)
d = 0:7;
totalP = (nWalls*(wallDensity + fluidDensity))+((nWalls/4)+1)^2;
tempPtCoord =[];
criticalP = points((1+(nWalls/8)*d),:);%%%hamishe bayad 8 ta bashe
direction = zeros(nWalls,3);
% tic
for y = 1:nWalls
    direction(y,:) = points(y,:) - center;
end
criticalDir = zeros(8,3);
criticalP1 = zeros(8,3);
criticalP2 = zeros(8,3);
for y = 1:8
    criticalDir(y,:) = (criticalP(y,:) - center);
    criticalP1(y,:) = center + criticalDir(y,:).*wall;
    if mod(y,2) == 1
        criticalP2(y,:) = center + criticalDir(y,:).*fluid;
    else
        criticalP2(y,:) = center + criticalDir(y,:).*fluid*(4.2/3.5);
    end
end
%%%lumen DOMAIN MESH
wallP = cell(1,wallDensity+1);
for b = wallDensity+1:-1:1
    wallP{b}(1:nWalls,:) = zeros(nWalls,3);
    for u = 1:nWalls
        temp = ((1-wall)/wallDensity)*(b-1);
        wallP{b}(u,:) = center + direction(u,:).*(wall+temp);
    end
end
%%%%line Point generation
counter = 0;
lineP = cell(1,4);
for h = 1:2:8
    counter = counter+1;
    DirLine = criticalP2(h,:)-center;
    lineP{counter}(1:(nWalls/8),:) = zeros((nWalls/8),3);
    for u = 1:(nWalls/8)
        lineP{counter}(u,:) = center + 8*DirLine*u/nWalls;
    end
end
% toc
cp2 =  [criticalP2;criticalP2(1,:)];
dir = zeros(8,3);
circle = cell(1,8);
square = cell(1,8);
dirSC = cell(1,8);
for b = 1:8
    circle{b} = zeros((nWalls/8),3);
    dir(b,:) = cp2(b+1,:) - cp2(b,:);
    temp = (1+(b-1)*(nWalls/8));
    circle{b} =  wallP{1}(temp:b*(nWalls/8),:);
    square{b} = zeros((nWalls/8),3);
    dirSC{b} = zeros((nWalls/8),3);
    for u = 1:(nWalls/8)
        temp = dir(b,:)*(8*(u-1)/nWalls);
        square{b}(u,:) = cp2(b,:) + temp;
        dirSC{b}(u,:) = circle{b}(u,:) - square{b}(u,:);
    end
end
fluidP = cell(1,8);
for h=1:8
    fluidP{h} = zeros(nWalls/8,3,(fluidDensity+1));
    for u = 1:nWalls/8
        for c = 1:(fluidDensity+1)
            temp = dirSC{h}(u,:)*((c-1)/fluidDensity);
            fluidP{h}(u,:,c) = square{h}(u,:) + temp;
            %                  scatter3(fluidP{h}(u,1,c),fluidP{h}(u,2,c),fluidP{h}(u,3,c))
        end
    end
end

inMesh = cell(1,wallDensity+fluidDensity+1);
squareP = cell(1,(nWalls/4)+1);
% toc
%%mini square1
temp = flipud(lineP{2});
terminal = [temp;center];%square{4};square{5}(1,:)];
temp = flipud(square{1});
start = [square{2}(1,:);temp];%flipud(lineP{2});center]
dir = terminal-start;
sPart1 = cell(1,(nWalls/8)+1);
sPart2 = cell(1,(nWalls/8)+1);
sPart4 = cell(1,(nWalls/8)+1);
sPart3 = cell(1,(nWalls/8)+1);
for k = 1:(nWalls/8)+1
    
    for u =1:(nWalls/8)+1
        sPart1{k}(u,:) = start(k,:)+ dir(k,:)*(u-1)/(nWalls/8);
        
    end
end
% toc
%%mini square2
terminal=[square{4};square{5}(1,:)];
start = [flipud(lineP{2});center];
dir = terminal-start;
for k = 1:(nWalls/8)+1
    
    for u = 1:(nWalls/8)
        sPart2{k}(u,:) = start(k,:) + dir(k,:)*(u)/(nWalls/8);
        
    end
end
% toc
%%mini square4
start = flipud(square{8});
terminal = lineP{4};
dir = terminal-start;
for k = 1:(nWalls/8)
    
    for u = 1:(nWalls/8)+1
        sPart4{k}(u,:) = start(k,:) + dir(k,:)*(u-1)/(nWalls/8);
        
    end
end
% toc
%%mini square3
start = lineP{4};
terminal = [square{5}(2:end,:);square{6}(1,:)];
dir = terminal-start;
for k = 1:(nWalls/8)
    for u =1:(nWalls/8)
        sPart3{k}(u,:) = start(k,:)+ dir(k,:)*(u)/(nWalls/8);
        
    end
end

for k=1:(nWalls/8)+1
    squareP{k} = [sPart1{k};sPart2{k}];
    % scatter3(squareP{k}(:,1),squareP{k}(:,2),squareP{k}(:,3))
end
for k= 1:(nWalls/8)
    squareP{k+((nWalls/8)+1)} = [sPart4{k};sPart3{k}];
    % scatter3(squareP{k+((nWalls/8)+1)}(:,1),squareP{k+((nWalls/8)+1)}(:,2),squareP{k+((nWalls/8)+1)}(:,3))
end

% toc
%%%csWritting

for k = 1:wallDensity
    inMesh{k+1} = wallP{(wallDensity+1)-k};
    
end

for k = (fluidDensity+1):-1:1
    f=[];
    for u = 1:8
        b = fluidP{u}(1:(nWalls/8),1:3,k);
        f = [f;b];
    end
    inMesh{fluidDensity+2+wallDensity-k} = f;
    
    
end
inMesh{1} = points;
indexMx=[];
% toc
for f = 1:(wallDensity+fluidDensity)
    inMesh{f} = [inMesh{f}(((nWalls/8)+1):end,:);inMesh{f}((1:(nWalls/8)),:)];
    ptIndex = ((ptCoordMxSize+1):(ptCoordMxSize+nWalls))';
    tempPtCoord = [tempPtCoord;inMesh{f}];
    inMesh{f} = [ptIndex,inMesh{f}(:,1:3)];
    indexMx =  [indexMx,ptIndex] ;
    ptCoordMxSize = ptCoordMxSize + nWalls;
end
for f = 1:(nWalls/4)+1
    ptIndex = ((ptCoordMxSize+1):(ptCoordMxSize+(nWalls/4)+1))';
    tempPtCoord = [tempPtCoord; squareP{f}];
    squareP{f} = [ptIndex,squareP{f}(:,1:3)];
    indexMx(1:(nWalls/4)+1,wallDensity+fluidDensity+f) = ptIndex ;
    ptCoordMxSize = ptCoordMxSize + (nWalls/4)+1;
end
ptCoordMx(ptCoordMxSize-totalP+1:(ptCoordMxSize),:) = tempPtCoord;




function [faceMx, bif,bifFaceMx,faceMxSize] = innerFaceMxMaker (~,faceMx, ~, bif,nWalls,~,~,~,~,bifFaceMx,csDensity,~,d,ptCoordMx,rightHandTest,ptCoordMxSize,faceMxSize)
n = faceMxSize.in;
offset = ptCoordMxSize-faceMxSize.in-d;
nS1 = faceMxSize.sur;
iFace1 = n+1;
nS1 = nS1+1;

%%%%%%%%%%%%faceMx developmet%%%%%%%%%%%%%%%%%%%
localFaceMxIn = zeros (d*2,6);
localFaceMxSur = zeros(nWalls,6);
iFace = 1;
nS = 1;
startFace = iFace;
startFace2 = nS;
for t = 1:csDensity-1
    for a = 1:(nWalls)
        
        %%%first & last wall
        if a == 1 && t ==1
            n = n+1;
            localFaceMxIn (iFace,:) = [n+offset n+offset+nWalls n+offset+d+nWalls n+d+offset n+offset n+offset+nWalls-1];
            localFaceMxSur(nS,:) = [n+offset n+d+offset n+1+d+offset n+1+offset n+offset 0];
            iFace = iFace + 1;
            nS = nS+1;
            %%%first step and last cell
        elseif t==1 && a == nWalls
            n = n+1;
            localFaceMxIn(iFace,:) = [n+offset n+offset+nWalls n+offset+nWalls+d n+offset+d n+offset n+offset-1];
            localFaceMxSur(nS,:) = [n+offset n+offset+d n+1-nWalls+offset+d n+1-nWalls+offset n+offset 0];
            iFace = iFace+1;
            nS = nS+1;
        elseif t == 1
            n = n+1;
            localFaceMxIn(iFace,:) = [n+offset n+offset+nWalls n+offset+d+nWalls n+d+offset n+offset n+offset-1];
            localFaceMxSur(nS,:) = [n+offset n+d+offset n+1+d+offset n+1+offset n+offset 0];
            iFace = iFace + 1;
            nS = nS+1;
        elseif a ==1
            n = n+1;
            localFaceMxIn(iFace,:) = [n+offset n+offset+nWalls n+offset+d+nWalls n+d+offset n+offset n+offset+nWalls-1];
            localFaceMxIn(iFace+1,:) = [n+offset n+d+offset n+1+d+offset n+1+offset n+offset n+offset-nWalls];
            iFace = iFace + 2;
        elseif a ==nWalls
            n = n+1;
            localFaceMxIn(iFace,:) = [n+offset n+offset+nWalls n+offset+nWalls+d n+offset+d n+offset n+offset-1];
            localFaceMxIn(iFace+1,:) = [n+offset n+offset+d n+1-nWalls+offset+d  n+1-nWalls+offset n+offset n+offset-nWalls];
            iFace = iFace+2;
        else
            n = n+1;
            localFaceMxIn(iFace,:) = [n+offset n+offset+nWalls n+offset+d+nWalls n+d+offset n+offset n+offset-1];
            localFaceMxIn(iFace+1,:) = [n+offset n+d+offset n+1+d+offset n+1+offset n+offset n+offset-nWalls];
            iFace = iFace + 2;
        end
    end
end

% betweeen square and circle part
start = n+offset+nWalls+1;
r = (nWalls/4);
a1 = start:(start+(r)-1);
b = (start+(r)):((r)+1):(start+((r)+1)*(r));
temp = (start+1+r*(r+1):start+r+r*(r+1))';
c = flipud(temp)';
tempN = (start+r+1:(r+1):(start+r*(r+1)))';
temp = flipud(tempN);
d1 = (temp)';
squareIndex = 1;
square = [a1 b c d1];
square = [square square(1)];
b2 = b-1;
c2 = c-((nWalls/4)+2);
d2 = d1-((nWalls/4)+1);
squCell = [a1 b2 c2 d2];
iSqu = 0 ;
for a = 1:(nWalls)
    
    %%%avalin step va alalin wall
    if a == 1
        n = n+1;
        iSqu = iSqu +1;
        localFaceMxIn(iFace,:) = [n+offset square(squareIndex) square(squareIndex)+d n+d+offset n+offset n+offset+nWalls-1];
        localFaceMxIn(iFace+1,:) = [n+offset n+d+offset n+1+d+offset n+1+offset n+offset n+offset-nWalls];
        localFaceMxIn(iFace+2,:) = [square(squareIndex) square(squareIndex+1) square(squareIndex+1)+d square(squareIndex)+d n+offset squCell(iSqu)];
        iFace = iFace+3;
        %%% avain step va akharin cell
    elseif a == nWalls
        n = n+1;
        iSqu = iSqu +1;
        squareIndex = squareIndex+1;
        localFaceMxIn(iFace,:) = [n+offset square(squareIndex) square(squareIndex)+d n+offset+d n+offset n+offset-1];
        localFaceMxIn(iFace+1,:) = [n+offset n+offset+d n+1-nWalls+offset+d  n+1-nWalls+offset n+offset n+offset-nWalls];
        localFaceMxIn(iFace+2,:) = [square(squareIndex) square(squareIndex+1) square(squareIndex+1)+d square(squareIndex)+d n+offset squCell(iSqu)];
        iFace = iFace+3;
    else
        squareIndex = squareIndex+1;
        iSqu = iSqu +1;
        n = n+1;
        localFaceMxIn(iFace,:) = [n+offset square(squareIndex) square(squareIndex)+d n+d+offset n+offset n+offset-1];
        localFaceMxIn(iFace+1,:) = [n+offset n+d+offset n+1+d+offset n+1+offset n+offset n+offset-nWalls];
        localFaceMxIn(iFace+2,:) = [square(squareIndex) square(squareIndex+1) square(squareIndex+1)+d square(squareIndex)+d n+offset squCell(iSqu)];
        iFace = iFace+3;
    end
end




for t = 1:(nWalls/4)-1
    
    for t2 = 1:(nWalls/4)-1
        n=n+1;
        localFaceMxIn(iFace,:) = [n+(nWalls/4)+1+offset n+(nWalls/4)+2+offset n+(nWalls/4)+2+offset+d n+(nWalls/4)+1+offset+d n+offset n+offset+(nWalls/4)+1];
        localFaceMxIn(iFace+1,:) = [n+(nWalls/4)+2+offset n+1+offset n+1+offset+d n+(nWalls/4)+2+offset+d n+offset n+offset+1];
        iFace = iFace+2;
    end
    n=n+1;
    localFaceMxIn(iFace,:) = [n+(nWalls/4)+1+offset n+(nWalls/4)+2+offset n+(nWalls/4)+2+offset+d n+(nWalls/4)+1+offset+d n+offset n+offset+(nWalls/4)+1];
    iFace = iFace+1;
    offset = offset +1;
end
for t=1:(nWalls/4)-1
    n = n+1;
    localFaceMxIn(iFace,:) = [n+(nWalls/4)+2+offset n+1+offset n+1+offset+d n+(nWalls/4)+2+offset+d n+offset n+offset+1];
    iFace = iFace+1;
end
endFace = iFace-1;
endFace2 = nS-1;
if rightHandTest == -1
    localFaceMxIn(startFace:endFace,1:6) = [localFaceMxIn(startFace:endFace,1:4),localFaceMxIn(startFace:endFace,6),localFaceMxIn(startFace:endFace,5)];
    localFaceMxSur(startFace2:endFace2,1:6) = [localFaceMxSur(startFace2:endFace2,1:4),localFaceMxSur(startFace2:endFace2,6),localFaceMxSur(startFace2:endFace2,5)];
end

localFaceMxIn = localFaceMxIn(any(localFaceMxIn,2),:);
tool = iFace1+length(localFaceMxIn)-1;
faceMx.in(iFace1:(tool),:) = localFaceMxIn;
faceMxSize.in = tool;
faceMx.sur(nS1:(nS1+nWalls-1),:) = localFaceMxSur;
faceMxSize.sur = (nS1+nWalls-1);


function [faceMx, bif,bifFaceMx,faceMxSize] = inletFaceMaker (~,faceMx, ~, bif,nWalls,~,~,~,~,bifFaceMx,csDensity,~,d,ptCoordMx,rightHandTest,j,ptCoordMxSize,faceMxSize)
% % bifIndex = 0;
if j == 1
    if  faceMxSize.inlet == 0
        n = 0;
        offset = ptCoordMxSize-d;
    else
        n = faceMxSize.inlet;
        offset =ptCoordMxSize-faceMxSize.inlet-d;
    end
end
if j ~=1
    n = faceMxSize.in;
    offset =ptCoordMxSize-faceMxSize.in-d;
end
iFace1 = n+1;
%%%%%%%%%%%%faceMx developmet%%%%%%%%%%%%%%%%%%%
localFaceMxInlet = zeros(64,6);
localFaceMxIn = zeros(d*2,6);
iFace = 1;
startFace = iFace;
if j == 1
    for t = 1:csDensity-1
        for a = 1:(nWalls)
            %%%avalin step va alalin wall
            if a == 1 && t ==1
                n = n+1;
                localFaceMxInlet(iFace,:) = [n+offset n+1+offset n+1+nWalls+offset n+nWalls+offset n+offset 0];
                iFace = iFace + 1;
                %%% avain step va akharin cell
            elseif t == 1 && a == nWalls
                n = n+1;
                localFaceMxInlet(iFace,:) = [n+offset  n+1+offset-nWalls n+1+offset n+offset+nWalls n+offset 0];
                iFace = iFace+1;
            elseif t == 1
                n = n+1;
                localFaceMxInlet(iFace,:) = [n+offset n+1+offset n+1+nWalls+offset n+nWalls+offset n+offset 0];
                iFace = iFace + 1;
            elseif a == 1
                n = n+1;
                localFaceMxInlet(iFace,:) = [n+offset n+1+offset n+1+nWalls+offset n+nWalls+offset n+offset 0];
                iFace = iFace + 1;
            elseif a == nWalls
                n = n+1;
                localFaceMxInlet(iFace,:) = [n+offset  n+1+offset-nWalls n+1+offset n+offset+nWalls n+offset 0];
                iFace = iFace+1;
            else
                n = n+1;
                localFaceMxInlet(iFace,:) = [n+offset n+1+offset n+1+nWalls+offset n+nWalls+offset n+offset 0];
                iFace = iFace + 1;
            end
        end
    end
end

if j~= 1
    for t = 1:csDensity-1
        for a = 1:(nWalls)
            
            %%%avalin step va alalin wall
            if a == 1 && t ==1
                n = n+1;
                localFaceMxIn(iFace,:) = [n+offset n+1+offset n+1+nWalls+offset n+nWalls+offset n+offset n+offset-d];
                iFace = iFace + 1;
                %%% avain step va akharin cell
            elseif t==1 && a == nWalls
                n = n+1;
                localFaceMxIn(iFace,:) = [n+offset  n+1+offset-nWalls n+1+offset n+offset+nWalls n+offset n+offset-d];
                iFace = iFace+1;
            elseif t == 1
                n = n+1;
                localFaceMxIn(iFace,:) = [n+offset n+1+offset n+1+nWalls+offset n+nWalls+offset n+offset n+offset-d];
                iFace = iFace + 1;
            elseif a ==1
                n = n+1;
                localFaceMxIn(iFace,:) = [n+offset n+1+offset n+1+nWalls+offset n+nWalls+offset n+offset n+offset-d];
                iFace = iFace + 1;
            elseif a ==nWalls
                n = n+1;
                localFaceMxIn(iFace,:) = [n+offset  n+1+offset-nWalls n+1+offset n+offset+nWalls n+offset n+offset-d];
                iFace = iFace+1;
            else
                n = n+1;
                localFaceMxIn(iFace,:) = [n+offset n+1+offset n+1+nWalls+offset n+nWalls+offset n+offset n+offset-d];
                iFace = iFace + 1;
            end
        end
    end
end
% betweeen square and circle part
start = n+offset+nWalls+1;
r = nWalls/4;
a1 = start:(start+r-1);
b = (start+r):(r+1):(start+(r+1)*r);
temp = ((start+1+r*(r+1):start+r+(r)*(r+1)))';
c = flipud(temp)';
temp = (start+r+1:(r+1):(start+(r)*(r+1)))';
d1 = (flipud(temp))';
squareIndex = 1;
square = [a1 b c d1];
square = [square square(1)];
iSqu = 0 ;
if j == 1
    for a = 1:(nWalls)
        
        %%%avalin step va alalin wall
        if a == 1
            n = n+1;
            iSqu = iSqu +1;
            localFaceMxInlet(iFace,:) = [n+offset n+1+offset square(squareIndex+1) square(squareIndex) n+offset 0];
            iFace = iFace+1;
            %%% avain step va akharin cell
        elseif a == nWalls
            n = n+1;
            iSqu = iSqu +1;
            squareIndex = squareIndex+1;
            localFaceMxInlet(iFace,:) = [n+offset  n+1+offset-nWalls square(squareIndex+1) square(squareIndex) n+offset 0];
            iFace = iFace+1;
        else
            squareIndex = squareIndex+1;
            iSqu = iSqu +1;
            n = n+1;
            localFaceMxInlet(iFace,:) = [n+offset n+1+offset square(squareIndex+1) square(squareIndex) n+offset 0];
            iFace = iFace+1;
        end
    end
end

if j ~=1
    for a = 1:(nWalls)
        
        %%%avalin step va alalin wall
        if a == 1
            n = n+1;
            iSqu = iSqu +1;
            localFaceMxIn(iFace,:) = [n+offset n+1+offset square(squareIndex+1) square(squareIndex) n+offset n+offset-d];
            iFace = iFace+1;
            %%% avain step va akharin cell
        elseif a == nWalls
            n = n+1;
            iSqu = iSqu +1;
            squareIndex = squareIndex+1;
            localFaceMxIn(iFace,:) = [n+offset  n+1+offset-nWalls square(squareIndex+1) square(squareIndex) n+offset n+offset-d];
            iFace = iFace+1;
        else
            squareIndex = squareIndex+1;
            iSqu = iSqu +1;
            n = n+1;
            localFaceMxIn(iFace,:) = [n+offset n+1+offset square(squareIndex+1) square(squareIndex) n+offset n+offset-d];
            iFace = iFace+1;
        end
    end
end
if j == 1
    for t = 1:(nWalls/4)-1
        
        for t2 = 1:(nWalls/4)-1
            n = n+1;
            localFaceMxInlet(iFace,:) = [n+offset n+1+offset n+(nWalls/4)+2+offset n+(nWalls/4)+1+offset n+offset 0];
            iFace = iFace+1;
        end
        n = n+1;
        localFaceMxInlet(iFace,:) = [n+offset n+1+offset n+(nWalls/4)+2+offset n+(nWalls/4)+1+offset n+offset 0];
        iFace = iFace+1;
        offset = offset +1;
    end
    for t=1:(nWalls/4)-1
        n = n+1;
        localFaceMxInlet(iFace,:) = [n+offset n+1+offset n+(nWalls/4)+2+offset n+(nWalls/4)+1+offset n+offset 0];
        iFace = iFace+1;
    end
    n = n+1;
    localFaceMxInlet(iFace,:) = [n+offset n+1+offset n+(nWalls/4)+2+offset n+(nWalls/4)+1+offset n+offset 0];
    endFace = iFace;
    if rightHandTest == -1
        localFaceMxInlet(startFace:endFace,1:6) = [localFaceMxInlet(startFace:endFace,1:4),localFaceMxInlet(startFace:endFace,6),localFaceMxInlet(startFace:endFace,5)];
    end
end
if j ~= 1
    for t = 1:(nWalls/4)-1
        
        for t2 = 1:(nWalls/4)-1
            n = n+1;
            localFaceMxIn(iFace,:) = [n+offset n+1+offset n+(nWalls/4)+2+offset n+(nWalls/4)+1+offset n+offset n+offset-d];
            iFace = iFace+1;
        end
        n = n+1;
        localFaceMxIn(iFace,:) = [n+offset n+1+offset n+(nWalls/4)+2+offset n+(nWalls/4)+1+offset n+offset n+offset-d];
        iFace = iFace+1;
        offset = offset +1;
    end
    for t=1:(nWalls/4)-1
        n = n+1;
        localFaceMxIn(iFace,:) = [n+offset n+1+offset n+(nWalls/4)+2+offset n+(nWalls/4)+1+offset n+offset n+offset-d];
        iFace = iFace+1;
    end
    n = n+1;
    localFaceMxIn(iFace,:) = [n+offset n+1+offset n+(nWalls/4)+2+offset n+(nWalls/4)+1+offset n+offset n+offset-d];
    endFace = iFace;
    if rightHandTest == -1
        localFaceMxIn(startFace:endFace,1:6) = [localFaceMxIn(startFace:endFace,1:4),localFaceMxIn(startFace:endFace,6),localFaceMxIn(startFace:endFace,5)];
    end
end
localFaceMxInlet = localFaceMxInlet(any(localFaceMxInlet,2),:);
localFaceMxIn = localFaceMxIn(any(localFaceMxIn,2),:);
tool = iFace1+length(localFaceMxInlet)-1;
tool1 = iFace1+length(localFaceMxIn)-1;
if tool ~= iFace1-1
    faceMx.inlet(faceMxSize.inlet+1:(tool),:) = localFaceMxInlet;
    faceMxSize.inlet = tool;
end
if tool1 ~=iFace1 -1
    faceMx.in(faceMxSize.in+1:(faceMxSize.in+length(localFaceMxIn)),:) = localFaceMxIn;
    faceMxSize.in = (faceMxSize.in+length(localFaceMxIn));
end

function rightHandTest = ifRightHand(point1,point2,center,velocity)
vector1 = point1-center;
vector2 = point2-center;
crossP = cross(vector1,vector2);
if dot(crossP,velocity) >= 0
    rightHandTest = 1;
else
    rightHandTest = -1;
end

function [separationPoints,ifMoreThan90] = computeSeparationPoints(c1,bifCoord,pointCoord,radius)
%each row in separationdirection is the normalized the S1, S2, S3 direction respectively.
Tangent = zeros(3);
radius = mean(radius,2);
for u=1:3
    Tangent(u,:)=normalize(c1(u,:)-bifCoord);
end
if Tangent(1,:) == 0;
    Tangent(1,:)=normalize(pointCoord(1,:)-bifCoord);
end
if Tangent(2,:) == 0;
    Tangent(2,:)=normalize(pointCoord(2,:)-bifCoord);
end
if Tangent(3,:) == 0;
    Tangent(3,:)=normalize(pointCoord(3,:)-bifCoord);
end
separationDirection1(1,:) = normalize((radius(2)*(Tangent(1,:)))+(radius(1)*(Tangent(2,:))));
separationDirection1(2,:) = normalize((radius(3)*(Tangent(1,:)))+(radius(1)*(Tangent(3,:))));
separationDirection1(3,:) = normalize((radius(3)*(Tangent(2,:)))+(radius(2)*(Tangent(3,:))));

separationDirection2(1,:) = normalize(((Tangent(1,:)))+((Tangent(2,:))));
separationDirection2(2,:) = normalize(((Tangent(1,:)))+((Tangent(3,:))));
separationDirection2(3,:) = normalize(((Tangent(2,:)))+((Tangent(3,:))));

separationDirection(1,:) = normalize(separationDirection1(1,:)+separationDirection2(1,:));
separationDirection(2,:) = normalize(separationDirection1(2,:)+separationDirection2(2,:));
separationDirection(3,:) = normalize(separationDirection1(3,:)+separationDirection2(3,:));
separationDirection = findOrientation(Tangent(1,:),Tangent(2,:),Tangent(3,:),c1(1,:),c1(2,:),c1(3,:),bifCoord,separationDirection);
% hold on
% axis equal
% scatter3 (c1(:,1),c1(:,2),c1(:,3),'filled','b')scatter3(bifCoord(:,1),bifCoord(:,2),bifCoord(:,3),'filled','r')
%%%find two most closest branches
angel(1) = computeAngle(c1(1,:),bifCoord,c1(2,:));
angel(2) = computeAngle(c1(2,:),bifCoord,c1(3,:));
angel(3) = computeAngle(c1(3,:),bifCoord,c1(1,:));
angleInDegrees = radtodeg(angel);
less = min(angel);
[~,b] = find(angel == less);
for i = 1:3
    %%%put dangerous check mark here
    if i ==3
        k = radius(1);
    else
        k = radius(i+1);
    end
    if angleInDegrees(i)<90
        m(i) = radius(i)/sin(atan(radius(i)/k));
        %         m(i) = radius(i)/sin(oneAngle);
        if angleInDegrees(i)<60 && radius(i)< k %%added 5/12/2017
            m(i) = (60/angleInDegrees(i))*(k/radius(i))*m(i);
        end
        ifMoreThan90(i) = 0;
    else
        m(i) = (radius(i)+k)/2;
        ifMoreThan90(i) = 1;
    end
    if m(i) > (radius(i)+k)
        m(i) = (radius(i)+k);
    end
end
separationPoints(1,:) = ((separationDirection(1,:))*(m(1)))+bifCoord;
separationPoints(2,:) = ((separationDirection(2,:))*(m(3)))+bifCoord;
separationPoints(3,:) = ((separationDirection(3,:))*(m(2)))+bifCoord;
% scatter3(separationPoints(1,1),separationPoints(1,2),separationPoints(1,3),'filled','r')% scatter3(separationPoints(2,1),separationPoints(2,2),separationPoints(2,3),'filled','g')% scatter3(separationPoints(3,1),separationPoints(3,2),separationPoints(3,3),'filled','b')
function separationDirection = findOrientation(mainVector,rightVector,leftVector,stemPoint,rightPoint,leftPoint,middlePoint,separationDirection)
mainVector = -mainVector;
cross1=cross((mainVector),(rightVector));
cross2=cross((mainVector),(leftVector));
if dot(cross1,cross2) >= 0
    %%compare angle and find the largest one
    angle1=computeAngle(rightPoint,middlePoint,stemPoint);
    angle2=computeAngle(leftPoint,middlePoint,stemPoint);
    if angle1 < angle2
        separationDirection(2,:) = -separationDirection(2,:);
    else
        separationDirection(1,:) = -separationDirection(1,:);
    end
    
end
%%%checing hte separation direction3 between the second and third branch,
%%%reorder main is converted to right,
vectors= [mainVector;rightVector;leftVector];
points = [stemPoint;rightPoint;leftPoint;middlePoint];
mainVector = -vectors(2,:); %%mainvector should always be flipped
rightVector = -vectors(1,:); %%because it has been fliped before the function
cross11 = cross((mainVector),(rightVector));
cross22 = cross((mainVector),(leftVector));
rightPoint = points(1,:);
stemPoint = points(2,:);
if dot(cross11,cross22) >= 0
    %%compare angle and find the largest one
    angle1=computeAngle(rightPoint,middlePoint,stemPoint);
    angle2=computeAngle(leftPoint,middlePoint,stemPoint);
    if angle1 < angle2
        separationDirection(3,:) = -separationDirection(3,:);
        
    else
        % separationDirection(2,:) = -separationDirection(2,:); %%I am not sure about it it may be doublicated
    end
    
end


function Angle = computeAngle(point_1,middle_point,point_2)
axis1=point_1-middle_point;
axis2=point_2-middle_point;
Angle=atan2(norm(cross(axis1,axis2)),dot(axis1,axis2));




function nSteps= computenSteps(arcLen,criteria,radius)
nSteps = zeros(3,1);
radius = mean(radius,2);
if isinteger(criteria)
    nSteps(1) = round((arcLen(1)/radius(1)))*(criteria);
    nSteps(2) = round((arcLen(2)/radius(2)))*(criteria);
    nSteps(3) = round((arcLen(3)/radius(3)))*(criteria);
else
    nSteps(1) = round(((arcLen(1)/radius(1)))*(criteria));
    nSteps(2) = round(((arcLen(2)/radius(2)))*(criteria));
    nSteps(3) = round(((arcLen(3)/radius(3)))*(criteria));
end
for i=1:3
    if nSteps(i) == 0 || nSteps(i) == 1 || nSteps(i) == 2
        nSteps(i) = 2;
    end
end
function wallDistribution = computeWallDistribution(controlPoints,separationPoints,bifCoord,nWalls)
%%the wallDistribution1 shows distribution of walls around bifucation
%%site,the distribution start from bifurcation axis in clock wise each
%%column of wallDistribution shows the wallDistribution from C1 in clock
%%wise for each vessels
halfWalls=nWalls/2;
criterion=pi/halfWalls;

wallDistribution(1)=round((computeAngle(separationPoints(1,:),bifCoord,controlPoints(1,:)))/criterion);
wallDistribution(2)=round((computeAngle(separationPoints(1,:),bifCoord,controlPoints(2,:)))/criterion);
wallDistribution(3)=round((computeAngle(separationPoints(2,:),bifCoord,controlPoints(2,:)))/criterion);
wallDistribution(4)=round((computeAngle(separationPoints(2,:),bifCoord,controlPoints(1,:)))/criterion);
wallDistribution(5)=round((computeAngle(separationPoints(3,:),bifCoord,controlPoints(2,:)))/criterion);
wallDistribution(6)=round((computeAngle(separationPoints(3,:),bifCoord,controlPoints(1,:)))/criterion);

if wallDistribution(1) == 0
    wallDistribution(1)= 1;
    wallDistribution(2)= wallDistribution(2)-1;
end

if wallDistribution(3) == 0
    wallDistribution(3)= 1;
    wallDistribution(4)= wallDistribution(4)-1;
end

if wallDistribution(6) == 0
    wallDistribution(6)= 1;
    wallDistribution(5)= wallDistribution(5)-1;
end

wallDistribution1=[wallDistribution(1);wallDistribution(2);wallDistribution(3);wallDistribution(4)];
wallDistribution2=[wallDistribution(1);wallDistribution(2);wallDistribution(5);wallDistribution(6)];
wallDistribution=cat(2,wallDistribution1,wallDistribution2,wallDistribution1);

function [Length] = computeArcLength(p0,p1,c0,c1, ~)
n=10;
t=0:(1/n):1;
bezier=zeros((n+1),3);
for i=1:(n+1)
    T=t(i);
    bezier(i,:)=((1-T)^3*p0+3*(T)*(1-T)^2*c0+3*(1-T)*(T)^2*c1+T^3*p1)';
end
Dist = zeros(n,1);
for j=1:n
    Dist(j) = norm(bezier(j+1,:)-bezier(j,:))';
end
Length=sum(Dist);
function bezierCurve = computeBezier(point1, c0, point2, c1,nSteps)
bezierCurve = zeros((nSteps+1),3);
t=0:(1/nSteps):1;
for i = 1:(nSteps+1)
    T = t(i);
    bezierCurve(i,:)=((1-T)^3*point1+3*(T)*(1-T)^2*c0+3*(1-T)*(T)^2*c1+T^3*point2)';
end




function startPoint2=computeStartPoint(point1,center1,tangent1,center2,tangent2,radius)
vector=point1-center1;
test=cross(vector,tangent1);
if test ~= 0
    
    startPoint2=ComputeProjectedPoint(point1,center1,center2,tangent2,radius);
else
    normalVector1=normalize(cross(vector,tangent1));
    normalVector2=normalize(cross(tangent2,normalVector1));
    startPoint2=center2+(normalVector2*radius);
end
function controlPoints = computeControlPoints(separationPoints,bifCoord,radius)
radius = mean(radius,2);
axis1_2=separationPoints(2,:)-separationPoints(1,:);
axis1_3=separationPoints(3,:)-separationPoints(1,:);
bifAxis(1,:)=normalize(cross(axis1_3,axis1_2));
bifAxis(2,:)=normalize(cross(axis1_2,axis1_3));

controlPoints(1,:)=bifCoord+bifAxis(1,:)*mean(radius);
controlPoints(2,:)=bifCoord+bifAxis(2,:)*mean(radius);

function bifurSitPoints=computeBifSitPointsType1(point1,point2,bifCoord,radius,radius1,radius2,wallDistribution)
counter= 0;
temp1 = point1-bifCoord;
tempAbs1 = (temp1(1)^2+temp1(2)^2+temp1(3)^2)^.5;

temp2 = point2-bifCoord;
tempAbs2 = (temp2(1)^2+temp2(2)^2+temp2(3)^2)^.5;
r = zeros (wallDistribution,1);
if tempAbs1 >= tempAbs2
    for t= 1:-(1/(wallDistribution)):0;
        counter= counter+1;
        r(counter,:)= [t 1-t]*[tempAbs1;tempAbs2];
    end
else
    for t= 1:-(1/(wallDistribution)):0;
        counter= counter+1;
        r(counter,:)= [t 1-t]*[tempAbs1;tempAbs2];
    end
end
counter1= 0;
b= normalize(point1-bifCoord);
b1= normalize(cross(b,(point2-bifCoord)));
b2= cross(b1,b);
bifurSitPoints = zeros(wallDistribution,3);
temp = computeAngle(point1,bifCoord,point2);
for h= 1:(wallDistribution);
    teta= (h-1)*(temp/wallDistribution);
    counter1= counter1+1;
    bifurSitPoints(h,1:3)=bifCoord + (cos(teta)*r(counter1)*b) + (r(counter1)*sin(teta)*b2);
    
end
function projectedPoint=ComputeProjectedPoint(point1,centerPoint1,centerPoint2,tangent,radius1)
direction1=normalize(point1-centerPoint1);
projectedVector=direction1-((dot(direction1,normalize(tangent)))*normalize(tangent));
projectedPoint=(normalize(projectedVector)*radius1)+centerPoint2;
function vector = normalize(aVector)
vector = aVector/norm(aVector);
