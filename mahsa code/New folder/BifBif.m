%% MATLAB code for Parametric Meshing of bifurcation_bifurcation
%% Mahsa 01/06/2014
%% Linninger research Group

function [ptCoordMx, faceMx, bif,faceMxSize,ptCoordMxSize] = BifBif(bzNWK,nwk,bif,ptCoordMx, faceMx,nWalls,criteria,faceMxSize,ptCoordMxSize,wallDensity,fluidDensity)
%  clc;
%  clear all;
%  close all;
csDensity = wallDensity + fluidDensity;
p = bzNWK.obj2msh(:,1) == 2;
Mx = bzNWK.obj2msh(p,(2:4));
for i = 1: size (Mx,1)
    
    tangent1=[];tangent2=[];connect1 = []; 
    spline = Mx(i,1);
    radius = bzNWK.diaSpline(spline,:);
    radius(radius ==0)=[];
    connectedSpline1 = Mx(i,2);connectedSpline2 = Mx(i,3);
    [r1, ~] = find(bif.circle(:,1) == connectedSpline1);
    circlePoints1 = bif.circle(r1,:);
    [r3, ~] = find(bif.index(:,1) == connectedSpline1);
    pIndex1 = (bif.index(r3:((wallDensity+fluidDensity+((nWalls)/4)+1)+r3-1),2:nWalls+1))';
    nBifBif = i;
    [r2,~] = find(bif.circle(:,1) == connectedSpline2);
    circlePoints2 = bif.circle(r2,:);
    [r4,~] = find(bif.index(:,1) == connectedSpline2);
    pIndex2 = (bif.index(r4:((csDensity+((nWalls)/4))+r4),2:nWalls+1))';
    circleP1 = zeros(nWalls,3);
    circleP2 = zeros(nWalls,3);
    for k = 1:(nWalls)
        f = 2+(k-1)*3;
        xx1 = circlePoints1(f);
        yy1 = circlePoints1(f+1);
        zz1 = circlePoints1(f+2);
        circleP1(k,:) = [xx1 yy1 zz1];
        xx2 = circlePoints2(f);
        yy2 = circlePoints2(f+1);
        zz2 = circlePoints2(f+2);
        circleP2(k,:) = [xx2 yy2 zz2];
    end
    allBz = bzNWK.bzSplines(spline,:);
    if allBz(1) == bzNWK.bzSplines(connectedSpline1,4)
        allBz = allBz(allBz~=0);
    else
        allBz = fliplr(allBz(allBz~=0));
    end
    nBz = length(allBz)/4;
    c0 = [];c1 = [];p0 = [];p1 = [];
    for j = 1:nBz
        p0 = [p0;bzNWK.bzPtCrdMx(allBz((j-1)*4 + 1),:)];
        c0 = [c0;bzNWK.bzPtCrdMx(allBz((j-1)*4 + 2),:)];
        c1 = [c1;bzNWK.bzPtCrdMx(allBz((j-1)*4 + 3),:)];
        p1 = [p1;bzNWK.bzPtCrdMx(allBz((j-1)*4 + 4),:)];
        tangent1(j,:) = normalize(c0(j,:) - p0(j,:));
        tangent2(j,:) = - normalize (c1(j,:) - p1(j,:));
    end
    
    %%%%%%%%%%%angular way
    
    tangent1(1,:) = normalize(cross(ptCoordMx(pIndex1(2,1),:)-p0(1,:),(ptCoordMx(pIndex1(1,1),:)-p0(1,:))));
    tangent2(end,:) = normalize(cross(ptCoordMx(pIndex2(2,1),:)-p1(end,:),(ptCoordMx(pIndex2(1,1),:)-p1(end,:))));
    test1 = normalize(c0(1,:)-p0(1,:));
    test2 = normalize(c1(end,:)-p1(end,:));
    
    % %  k1=-1 means the counter-clock rotation of point in
    % input when looking from the bif-bif part
    if dot(tangent1(1,:),test1) <= 0
        k1 = -1;
        tangent1(2:(end),:) = -tangent1(2:(end),:);
        tangent2(1:(end-1),:) = -tangent2(1:(end-1),:);
    else
        k1 = 1;
    end
    % % %k2=-1 means the coutner-clock rotaion of point in output when looking from outside of bif-bif part
    if dot(tangent2(end,:),test2) >= 0
        k2 = -1;
    else
        k2 = 1;
    end
    if (dot(tangent1(1,:),test1) >= 0 && dot(tangent2(end,:),test2) >= 0) || (dot(tangent1(1,:),test1) <= 0 && dot(tangent2(end,:),test2) <= 0)
        tangent2(end,:) = -tangent2(end,:);
        k = -1;
    else
        k = 1;
    end
    
    %%%%

    if length(radius) ~= nBz+1
        radius(end:nBz+1) = radius(end);
    end

    first = ptCoordMx(pIndex1(1,1),:);
    projectP = zeros(nBz,3);
    for j = 1:nBz
        projectP(j,:) = computeStartPoint(first,p0(j,:),tangent1(j,:),p1(j,:),tangent2(j,:),radius(j+1));
        first = projectP(j,:);
    end
    
    %%%%%compute starting point
    
    group2 = pIndex2(1,1):nWalls/4:pIndex2(end,1);
    distance = zeros(1,length(group2));
    for n = 1:length(group2);
        distance(n) = computeDistance(projectP(end,:),ptCoordMx(group2(n),:));
    end
    minDist = min(distance);
    [~,cols] = find(distance == minDist);
    matchP = group2(cols);
    [row,~] = find(pIndex2 == matchP);
    cirNew = [];
    for i = 1:(csDensity)
        cirNew = [cirNew,[pIndex2(row:end,i);pIndex2(1:row,i)]];
    end
    
    
    if k == -1
        cirNew = flipud(cirNew);
    end
    Err = computeAngle(projectP(end,:),p1(end,:),ptCoordMx(matchP,:));
    
    if Err > 1.57075
        Err = 3.1415 - Err;
    end
    square = pIndex2(1:(nWalls/4)+1,csDensity+1:csDensity+(nWalls/4)+1);
    squareCell = square(1:nWalls/4,1:end-1);
    
    if row == (nWalls/4)+1
        if k == -1
            square = flipud(square);
            squareCell = flipud(squareCell);
        else
            square = flipud(square)';
            squareCell = flipud(squareCell)';
        end
    elseif row == (nWalls/2)+1
        if k == -1
            square = rot90(square,2)';
            squareCell = rot90(squareCell,2)';
        else
            square = rot90(square,2);
            squareCell = rot90(squareCell,2);
        end
    elseif row == 3*(nWalls/4)+1
        if k == -1
            square = fliplr(square);
            squareCell = fliplr(squareCell);
        else
            square = fliplr(square)';
            squareCell = fliplr(squareCell)';
        end
    elseif (row == 1 && k ==-1)
        square = square';
        squareCell = squareCell';
    end

    if k == -1
        cellIndex = cirNew;
    else
        cellIndex = [cirNew(nWalls,1:csDensity);cirNew(1:nWalls,1:csDensity)];
    end
    cirNew(1:(nWalls/4)+1,csDensity+1:csDensity+(nWalls/4)+1)=square;
    pIndex2 = cirNew;
    cellIndex(1:(nWalls/4),csDensity+1:csDensity+(nWalls/4)) = squareCell;
    circleP1 = ptCoordMx(pIndex1(:,1),:);
    
    
    if dot(normalize(projectP(end,:)-ptCoordMx(matchP,:)), normalize(ptCoordMx(pIndex2(2,1),:)-ptCoordMx(pIndex2(1,1),:))) <0
        Err = -1*Err;
    end
    radii=radius;
    arcLen = zeros(1,nBz);
    nSteps = zeros(1,nBz);
    for m=1:nBz
        arcLen(m) = computeArcLength(p0(m,:),p1(m,:),c0(m,:),c1(m,:));
        nSteps(m) = computenSteps(arcLen(m),criteria,mean([radius(m) radius(m+1)]));
    end
    error = Err/sum(nSteps);
    for m=1:nBz
        if  m ==1
            downP = pIndex1;
        end
        [ptCoordMx, faceMx, firstPoints,connect1, downP,k,faceMxSize,ptCoordMxSize] = computePtCoordMx([p0(m,:);p1(m,:)],c0(m,:),c1(m,:),nWalls,nSteps(m),radii(m),radii(m+1),circleP1,error,connect1,ptCoordMx, faceMx, bif, m, nBz,connectedSpline1,connectedSpline2,cols,pIndex1,pIndex2,wallDensity,fluidDensity, downP,cellIndex,k1,nBifBif,k2,k,faceMxSize,ptCoordMxSize,arcLen(m));
        circleP1 = firstPoints;
        %%%calculate new error
        initialP = circleP1((nWalls/8)+1,:);
        if m ~=nBz
            projectP1 = computeStartPoint(initialP,p0(m+1,:),tangent2(m,:),p1(end,:),tangent2(end,:),radius(end));
            newError = computeAngle(projectP1,p1(end,:),ptCoordMx(matchP,:));
            if dot((projectP1-ptCoordMx(matchP,:)), (ptCoordMx(pIndex2(2,1),:)-ptCoordMx(pIndex2(1,1),:))) <0
                newError = -1*newError;
            end
            error = newError/(sum(nSteps)-(sum(nSteps(1:m))));
        end
    end
end





function [ptCoordMx, faceMx, firstPoints,connect1, downP, k,faceMxSize,ptCoordMxSize] = computePtCoordMx(pointCoord,c0,c1,nWalls,nSteps,radius1,radius2,circleP,error,connect1,ptCoordMx, faceMx, ~,m, nBz,~,~,~,pIndex1,pIndex2,wallDensity,fluidDensity, downP, cellIndex,k1,nBifBif,k2,k,faceMxSize,ptCoordMxSize,arcLen)
position = computeBezier(pointCoord(1,:), c0, pointCoord(2,:), c1,nSteps,arcLen);
velocity = normalize(computeVelocity(pointCoord(1,:), c0, pointCoord(2,:), c1,nSteps));

if m == 1
    
    if dot((c0(1,:)-pointCoord(1,:)),normalize(cross((circleP(2,:)-pointCoord(1,:)),circleP(1,:)-pointCoord(1,:))))>=0
        k = -1;
    else
        
        k = 1;
    end
end
counter = 0;
r = zeros(1,nSteps+1);
for t=0:(1/(nSteps)):1;
    counter=counter+1;
    r(counter)=[1-t t]*[radius1;radius2];
end
if m==1
    circleP = [circleP((end-((nWalls/8)-1)):end,:) ; circleP(1:(end-(nWalls/8)),:)];
    wallPoints(:,1:3,1) = circleP;
else
    wallPoints(:,1:3,1) =   circleP;
end
%       scatter3(circleP(1,1),circleP(1,2),circleP(1,3),'filled','b');
%       scatter3(circleP(2,1),circleP(2,2),circleP(2,3),'filled','g');
projectedPoints = zeros(nSteps+1,3);
for i = 1:(nSteps+1)
    projectedPoints(i,:) = computeStartPoint(circleP(1,:),pointCoord(1,:),velocity(1,:),position(i,:),velocity(i,:),r(i));
    %      scatter3(projectedPoints(i,1),projectedPoints(i,2),projectedPoints(i,3));
end

for i = 2:(nSteps+1)
    b = normalize(projectedPoints(i,:)-position(i,:));
    counter = 0;
    b1 = normalize(velocity(i,:));
    b2 = cross(b1,b);
    
    for h = 1:nWalls
        teta = k*((h-1)*(2*pi/nWalls))+((-k)*(((i-1)*error)));
        counter = counter+1;
        wallPoints(h,1:3,i) = position(i,:)+(cos(teta)*r(i)*b)+(r(i)*sin(teta)*b2);
    end
    %   wallPoints(1:nWalls,1:3,i) =  [wallPoints((end-1):end,1:3,i); wallPoints(1:(end-2),1:3,i)];
%     plot3([wallPoints(:,1,i);wallPoints(1,1,i)],[wallPoints(:,2,i);wallPoints(1,2,i)],[wallPoints(:,3,i);wallPoints(1,3,i)],'linewidth',1,'color','g')
end
firstPoints = wallPoints(:,1:3,(nSteps+1));
%innermesh = computeInnerGrid(wallPoints(:,:,i),position(i,:),b1,teta,r(i), nWalls)
%%%ptCoordMx AND FaceMx development
wall = 0.8;
fluid = 0.45;
csDensity = wallDensity+fluidDensity;
d = (nWalls*(csDensity))+((nWalls/4)+1)^2;
for j = 1:nSteps
    xx = reshape(wallPoints(1:nWalls,1,j+1),1,nWalls);
    yy = reshape(wallPoints(1:nWalls,2,j+1),1,nWalls);
    zz = reshape(wallPoints(1:nWalls,3,j+1),1,nWalls);
    %                plot3(xx,yy,zz,'k')
    if j == nSteps && m ==nBz
        %              [inMesh, squareP,ptCoordMx,upP] = computeInnermesh([xx' yy' zz'],position(j+1,:),wallDensity,fluidDensity,wall,fluid,nWalls,ptCoordMx);
        
        if m ~= nBz
            [inMesh, squareP,ptCoordMx,upP] = computeInnermesh([xx' yy' zz'],position(j+1,:),wallDensity,fluidDensity,wall,fluid,nWalls,ptCoordMx);
            [faceMx] = innerFaceMxMaker (inMesh,faceMx,nWalls,csDensity,squareP,d,ptCoordMx,pIndex1,upP,m,downP,j,k1,nBifBif);
        else
            upP = pIndex2;
            %                         for t = upP:upP+d-1
            %                         scatter3(ptCoordMx(t,1),ptCoordMx(t,2),ptCoordMx(t,3))
            %                         end
            %                         [upP] = matchBif(downP, upP
            %                      [faceMx] = outletFaceMaker(faceMx,nWalls,csDensity,d,ptCoordMx,downP,upP);
            [faceMx,  faceMxSize] = matchFacebif2bifMaker (faceMx,nWalls,csDensity,d,ptCoordMx,pIndex1,upP,m,downP,j,cellIndex,k1,nBifBif,k2,  faceMxSize, ptCoordMxSize);
            %                      for t = downP:downP+d-1
            %     scatter3(ptCoordMx(t,1),ptCoordMx(t,2),ptCoordMx(t,3),'*k')
            % end
        end
        downP = upP;
    else
        [inMesh, squareP,ptCoordMx,upP, ptCoordMxSize] = computeInnermesh([xx' yy' zz'],position(j+1,:),wallDensity,fluidDensity,wall,fluid,nWalls,ptCoordMx, ptCoordMxSize);
        % for t = downP:downP+d-1
        %     scatter3(ptCoordMx(t,1),ptCoordMx(t,2),ptCoordMx(t,3))
        % end
        % for t = upP:upP+d-1
        %     scatter3(ptCoordMx(t,1),ptCoordMx(t,2),ptCoordMx(t,3))
        % end
        if m == 1 && j ==1
            [faceMx,  faceMxSize] = innerAttachFaceMxMaker (inMesh,faceMx,nWalls,csDensity,squareP,d,ptCoordMx,pIndex1,upP,k1,nBifBif,  faceMxSize, ptCoordMxSize);
            downP = upP;
        else
            [faceMx,  faceMxSize] = innerFaceMxMaker (inMesh,faceMx,nWalls,csDensity,squareP,d,ptCoordMx,pIndex1,upP,m,downP,j,k1,nBifBif,  faceMxSize, ptCoordMxSize);
            downP = upP;
        end
    end
end
function [faceMx,  faceMxSize] = innerAttachFaceMxMaker (~,faceMx,nWalls,csDensity,~,d,~,pIndex,upP,k1,nBifBif,  faceMxSize, ptCoordMxSize)
n = faceMxSize.in;
nS1 = faceMxSize.sur;
iFace1 = n+1;
nS1 = nS1+1;


localFaceMxIn = zeros (d*2,6);
localFaceMxSur = zeros(nWalls,6);
iFace = 1;
nS = 1;
startFace = iFace;
startFace2 = nS;
downP = pIndex;
offset = d*(nBifBif-1);
downP(end+1,1:csDensity)=downP(1,1:csDensity);
upP(end+1,1:csDensity)=upP(1,1:csDensity);
CSupP(1,:) =[upP(1:(nWalls/4),csDensity+1)',upP(((nWalls)/4)+1,csDensity+1:end),flipud(upP(1:(nWalls/4),end))',fliplr(upP(1,csDensity+1:end-1))];
CSupP(2,:) =[upP(1:(nWalls/4),csDensity+1)',upP(((nWalls)/4),csDensity+1:end-1),flipud(upP(1:(nWalls/4),end-1))',fliplr(upP(1,csDensity+1:end-1)),0];
%%%%%%%%%%%%faceMx developmet%%%%%%%%%%%%%%%%%%%
CSdownP(1,:) = [pIndex(1:(nWalls/4),csDensity+1)',pIndex(((nWalls)/4)+1,csDensity+1:end),flipud(pIndex(1:(nWalls/4),end))',fliplr(pIndex(1,csDensity+1:end-1))];
CSdownP(2,:) = [pIndex(1: (nWalls/4),csDensity+1)',pIndex(((nWalls)/4),csDensity+1:end-1),flipud(pIndex(1:(nWalls/4),end-1))',fliplr(pIndex(1,csDensity+1:end-1)),0];
for t = 1:csDensity
    i = 0;
    for a = 1:(nWalls)
        %%%alalin wall
        if t ==1
            i = i+1;
            localFaceMxIn(iFace,:) = [downP(i,t) downP(i+1,t)  downP(i+1,t+1) downP(i,t+1) downP(i,t) upP(i,t)+offset];
            localFaceMxIn(iFace+1,:) = [upP(i,t) upP(i+1,t) upP(i+1,t+1) upP(i,t+1) upP(i,t)+offset upP(i,t)+d+offset];
            localFaceMxSur(nS,:) = [upP(i,t) upP(i+1,t) downP(i+1,t) downP(i,t) 0 upP(i,t)+offset ];
            iFace = iFace + 2;
            nS = nS+1;
            if a ==1
                localFaceMxIn(iFace,:) = [upP(i,t) upP(i,t+1) downP(i,t+1) downP(i,t) upP(i,t)+offset upP(i,t)+nWalls-1+offset];
            else
                localFaceMxIn(iFace,:) = [upP(i,t) upP(i,t+1) downP(i,t+1) downP(i,t) upP(i,t)+offset upP(i,t)-1+offset];
            end
            iFace = iFace + 1;
        elseif t == csDensity
            i=i+1;
            localFaceMxIn(iFace,:) = [downP(i,t) downP(i+1,t) CSdownP(1,i+1) CSdownP(1,i) downP(i,t) upP(i,t)+offset];
            localFaceMxIn(iFace+1,:) =  [upP(i,t) upP(i+1,t) CSupP(1,i+1) CSupP(1,i) upP(i,t)+offset upP(i,t)+d+offset];
            if a == 1
                localFaceMxIn(iFace+2,:) =  [upP(i,t) CSupP(1,i) CSdownP(1,i) downP(i,t) upP(i,t)+offset upP(i,t)+nWalls-1+offset];
            else
                localFaceMxIn(iFace+2,:) =  [upP(i,t) CSupP(1,i) CSdownP(1,i) downP(i,t) upP(i,t)+offset upP(i-1,t)+offset];
            end
            localFaceMxIn(iFace+3,:) =  [upP(i,t) upP(i+1,t) downP(i+1,t) downP(i,t) upP(i,t-1)+offset upP(i,t)+offset];
            localFaceMxIn(iFace+4,:) =  [CSupP(1,i) CSupP(1,i+1) CSdownP(1,i+1) CSdownP(1,i) upP(i,t)+offset CSupP(2,i)+offset];
            iFace = iFace + 5;
        elseif (t > 1 && t < csDensity)
            i = i+1;
            localFaceMxIn(iFace,:) = [downP(i,t) downP(i+1,t) downP(i+1,t+1) downP(i,t+1) downP(i,t) upP(i,t)+offset];
            localFaceMxIn(iFace+1,:) = [upP(i,t) upP(i+1,t) upP(i+1,t+1) upP(i,t+1) upP(i,t)+offset upP(i,t)+d+offset];
            localFaceMxIn(iFace+2,:) = [upP(i,t) upP(i+1,t) downP(i+1,t) downP(i,t) upP(i,t-1)+offset upP(i,t)+offset ];
            if a ==1
                localFaceMxIn(iFace+3,:) = [upP(i,t) upP(i,t+1) downP(i,t+1) downP(i,t) upP(i,t)+offset upP(i,t)+nWalls-1+offset];
            else
                localFaceMxIn(iFace+3,:) = [upP(i,t) upP(i,t+1) downP(i,t+1) downP(i,t) upP(i,t)+offset upP(i,t)-1+offset];
            end
            iFace = iFace + 4;
        end
    end
end
for t = csDensity+1:csDensity+(nWalls/4)
    i=0;
    if  t ~= csDensity+(nWalls/4)
        for a = 1:(nWalls/4)
            if a == (nWalls/4)
                i=i+1;
                localFaceMxIn(iFace,:) = [downP(i,t) downP(i+1,t) downP(i+1,t+1) downP(i,t+1) downP(i,t) upP(i,t)+offset];
                localFaceMxIn(iFace+1,:) = [upP(i,t) upP(i+1,t) upP(i+1,t+1) upP(i,t+1) upP(i,t)+offset upP(i,t)+d+offset];
                localFaceMxIn(iFace+2,:) = [upP(i,t+1) upP(i+1,t+1) downP(i+1,t+1) downP(i,t+1) upP(i,t)+offset upP(i,t+1)+offset];
                iFace = iFace+3;
            else
                i=i+1;
                localFaceMxIn(iFace,:) = [downP(i,t) downP(i+1,t) downP(i+1,t+1) downP(i,t+1) downP(i,t) upP(i,t)+offset];
                localFaceMxIn(iFace+1,:) = [upP(i,t) upP(i+1,t) upP(i+1,t+1) upP(i,t+1) upP(i,t)+offset upP(i,t)+d+offset];
                localFaceMxIn(iFace+2,:) = [upP(i+1,t) upP(i+1,t+1) downP(i+1,t+1) downP(i+1,t) upP(i+1,t)+offset upP(i,t)+offset];
                localFaceMxIn(iFace+3,:) = [upP(i,t+1) upP(i+1,t+1) downP(i+1,t+1) downP(i,t+1) upP(i,t)+offset upP(i,t+1)+offset];
                iFace = iFace+4;
            end
        end
    else
        for a = 1:(nWalls/4)
            if a == (nWalls/4)
                i=i+1;
                localFaceMxIn(iFace,:) = [downP(i,t) downP(i+1,t) downP(i+1,t+1) downP(i,t+1) downP(i,t) upP(i,t)+offset];
                localFaceMxIn(iFace+1,:) = [upP(i,t) upP(i+1,t) upP(i+1,t+1) upP(i,t+1) upP(i,t)+offset upP(i,t)+d+offset];
                iFace = iFace+2;
            else
                i=i+1;
                localFaceMxIn(iFace,:) = [downP(i,t) downP(i+1,t) downP(i+1,t+1) downP(i,t+1) downP(i,t) upP(i,t)+offset];
                localFaceMxIn(iFace+1,:) = [upP(i,t) upP(i+1,t) upP(i+1,t+1) upP(i,t+1) upP(i,t)+offset upP(i,t)+d+offset];
                localFaceMxIn(iFace+2,:) = [upP(i+1,t) upP(i+1,t+1) downP(i+1,t+1) downP(i+1,t) upP(i+1,t)+offset upP(i,t)+offset];
                iFace = iFace+3;
                
                
            end
        end
        
    end
    
end
endFace = iFace-1;
endFace2 = nS-1;
if k1 == -1
    localFaceMxIn(startFace:endFace,1:6) = [localFaceMxIn(startFace:endFace,1:4),localFaceMxIn(startFace:endFace,6),localFaceMxIn(startFace:endFace,5)];
    localFaceMxSur(startFace2:endFace2,1:6) = [localFaceMxSur(startFace2:endFace2,1:4),localFaceMxSur(startFace2:endFace2,6),localFaceMxSur(startFace2:endFace2,5)];
end
localFaceMxIn = localFaceMxIn(any(localFaceMxIn,2),:);
tool = iFace1+length(localFaceMxIn)-1;
faceMx.in(iFace1:(tool),:) = localFaceMxIn;
faceMxSize.in = tool;
faceMx.sur(nS1:(nS1+nWalls-1),:) = localFaceMxSur;
faceMxSize.sur = (nS1+nWalls-1);
function [faceMx,  faceMxSize] = matchFacebif2bifMaker (faceMx,nWalls,csDensity,d,ptCoordMx,pIndex,upP,m,downP,j,cellIndex,k1,nBifBif,k2,  faceMxSize, ptCoordMxSize)
n = faceMxSize.in;
offset = ptCoordMxSize-faceMxSize.in-d;
nS1 = faceMxSize.sur;

iFace1 = n+1;
nS1 = nS1+1;
localFaceMxIn = zeros (d*2,6);
localFaceMxSur = zeros(nWalls,6);
iFace = 1;
nS = 1;
startFace = iFace;
startFace2 = nS;
upC = cellIndex;
offset = (nBifBif-1)*d;
upP(end+1,1:csDensity)=upP(1,1:csDensity);
CSupP(1,:) =[upP(1:(nWalls/4),csDensity+1)',upP(((nWalls)/4)+1,csDensity+1:end),flipud(upP(1:(nWalls/4),end))',fliplr(upP(1,csDensity+1:end-1))];
CSupP(2,:) =[upP(1:(nWalls/4),csDensity+1)',upP(((nWalls)/4),csDensity+1:end-1),flipud(upP(1:(nWalls/4),end-1))',fliplr(upP(1,csDensity+1:end-1)),0];
%%%%%%%%%%%%faceMx developmet%%%%%%%%%%%%%%%%%%%
if m == 1 && j == 1
    downP = pIndex;
    CSdownP(1,:) = [pIndex(1:(nWalls/4),csDensity+1)',pIndex(((nWalls)/4)+1,csDensity+1:end),flipud(pIndex(1:(nWalls/4),end))',fliplr(pIndex(1,csDensity+1:end-1))];
    CSdownP(2,:) = [pIndex(1: (nWalls/4),csDensity+1)',pIndex(((nWalls)/4),csDensity+1:end-1),flipud(pIndex(1:(nWalls/4),end-1))',fliplr(pIndex(1,csDensity+1:end-1)),0];
else
    CSdownP(1,:) = [downP(1:(nWalls/4),csDensity+1)',downP(((nWalls)/4)+1,csDensity+1:end),flipud(downP(1:(nWalls/4),end))',fliplr(downP(1,csDensity+1:end-1))];
    CSdownP(2,:) = [downP(1: (nWalls/4),csDensity+1)',downP(((nWalls)/4),csDensity+1:end-1),flipud(downP(1:(nWalls/4),end-1))',fliplr(downP(1,csDensity+1:end-1)),0];
end
downC = [downP(end,1:csDensity);downP(1:nWalls,1:csDensity)];
downC(1:nWalls/4,csDensity+1:csDensity+nWalls/4)=downP(1:nWalls/4,csDensity+1:csDensity+nWalls/4);
downC = downC+((nBifBif-1)*d);
middleC = downC+d;
downP(end+1,1:csDensity) = downP(1,1:csDensity);
%     for t=1:4
%         scatter3(ptCoordMx(upP(t,1),1),ptCoordMx(upP(t,1),2),ptCoordMx(upP(t,1),3),'*k')
%         scatter3(ptCoordMx(downP(t,1),1),ptCoordMx(downP(t,1),2),ptCoordMx(downP(t,1),3),'*r')
%     end
%     for t=1:3
%         scatter3(ptCoordMx(upP(t,4),1),ptCoordMx(upP(t,4),2),ptCoordMx(upP(t,4),3),'*k')
%         scatter3(ptCoordMx(downP(t,4),1),ptCoordMx(downP(t,4),2),ptCoordMx(downP(t,4),3),'*r')
%     end
for t = 1:csDensity
    i = 0;
    
    iCell = 1;
    
    for a = 1:(nWalls)
        if t ==1
            i = i+1;
            iCell = iCell+1;
            localFaceMxIn(iFace,:) = [upP(i,t) upP(i+1,t) upP(i+1,t+1) upP(i,t+1) middleC(iCell,t) upC(iCell,t)];
            localFaceMxIn(iFace+1,:) = [downP(i,t) downP(i+1,t)  downP(i+1,t+1) downP(i,t+1) downC(iCell,t) middleC(iCell,t)];
            localFaceMxIn(iFace+2,:) = [upP(i,t) upP(i,t+1) downP(i,t+1) downP(i,t) middleC(iCell,t) middleC(iCell-1,t)];
            localFaceMxSur(nS,:) = [upP(i,t) upP(i+1,t) downP(i+1,t) downP(i,t) 0 middleC(iCell,t) ];
            iFace = iFace + 3;
            nS = nS+1;
        elseif t == csDensity
            i=i+1;
            iCell = iCell+1;
            localFaceMxIn(iFace,:) =  [upP(i,t) upP(i+1,t) CSupP(1,i+1) CSupP(1,i) middleC(iCell,t) upC(iCell,t)];
            localFaceMxIn(iFace+1,:) = [downP(i,t) downP(i+1,t) CSdownP(1,i+1) CSdownP(1,i) downC(iCell,t) middleC(iCell,t)];
            localFaceMxIn(iFace+2,:) =  [upP(i,t) CSupP(1,i) CSdownP(1,i) downP(i,t) middleC(iCell,t) middleC(iCell-1,t)];
            localFaceMxIn(iFace+3,:) =  [upP(i,t) upP(i+1,t) downP(i+1,t) downP(i,t) middleC(iCell,t-1) middleC(iCell,t)];
            localFaceMxIn(iFace+4,:) =  [CSupP(1,i) CSupP(1,i+1) CSdownP(1,i+1) CSdownP(1,i) middleC(iCell,t) CSdownP(2,i)+d+offset];
            iFace = iFace + 5;
        elseif (t > 1 && t < csDensity)
            i = i+1;
            iCell = iCell+1;
            localFaceMxIn(iFace,:) = [upP(i,t) upP(i+1,t) upP(i+1,t+1) upP(i,t+1) middleC(iCell,t) upC(iCell,t)];
            localFaceMxIn(iFace+1,:) = [downP(i,t) downP(i+1,t)  downP(i+1,t+1) downP(i,t+1) downC(iCell,t) middleC(iCell,t)];
            localFaceMxIn(iFace+2,:) = [upP(i,t) upP(i+1,t) downP(i+1,t) downP(i,t) middleC(iCell,t-1) middleC(iCell,t) ];%%check
            localFaceMxIn(iFace+3,:) = [upP(i,t) upP(i,t+1) downP(i,t+1) downP(i,t) middleC(iCell,t) middleC(iCell-1,t)];
            iFace = iFace + 4;
        end
    end
end
for t = csDensity+1:csDensity+(nWalls/4)
    i = 0;
    if  t ~= csDensity+(nWalls/4)
        for a = 1:(nWalls/4)
            if a == (nWalls/4)
                i=i+1;
                localFaceMxIn(iFace,:) = [upP(i,t) upP(i+1,t) upP(i+1,t+1) upP(i,t+1)  middleC(i,t) upC(i,t)];
                localFaceMxIn(iFace+1,:) = [upP(i,t+1) upP(i+1,t+1) downP(i+1,t+1) downP(i,t+1) middleC(i,t) middleC(i,t+1)]; %%%%
                iFace = iFace+2;
            else
                i=i+1;
                localFaceMxIn(iFace,:) = [upP(i,t) upP(i+1,t) upP(i+1,t+1) upP(i,t+1)  middleC(i,t) upC(i,t)];
                localFaceMxIn(iFace+1,:) = [upP(i+1,t) upP(i+1,t+1) downP(i+1,t+1) downP(i+1,t) middleC(i+1,t) middleC(i,t)];
                localFaceMxIn(iFace+2,:) = [upP(i,t+1) upP(i+1,t+1) downP(i+1,t+1) downP(i,t+1) middleC(i,t) middleC(i,t+1)];
                iFace = iFace+3;
            end
        end
    else
        for a = 1:(nWalls/4)
            if a == (nWalls/4)
                i=i+1;
                localFaceMxIn(iFace,:) = [upP(i,t) upP(i+1,t) upP(i+1,t+1) upP(i,t+1)  middleC(i,t) upC(i,t)];
                iFace = iFace+1;
            else
                i=i+1;
                localFaceMxIn(iFace,:) = [upP(i,t) upP(i+1,t) upP(i+1,t+1) upP(i,t+1)  middleC(i,t) upC(i,t)];
                localFaceMxIn(iFace+1,:) = [upP(i+1,t) upP(i+1,t+1) downP(i+1,t+1) downP(i+1,t) middleC(i+1,t) middleC(i,t)];
                iFace = iFace+2;
            end
        end
    end
end
endFace = iFace-1;
endFace2 = nS-1;
if k1 == -1
    localFaceMxIn(startFace:endFace,1:6) = [localFaceMxIn(startFace:endFace,1:4),localFaceMxIn(startFace:endFace,6),localFaceMxIn(startFace:endFace,5)];
    localFaceMxSur(startFace2:endFace2,1:6) = [localFaceMxSur(startFace2:endFace2,1:4),localFaceMxSur(startFace2:endFace2,6),localFaceMxSur(startFace2:endFace2,5)];
end
localFaceMxIn = localFaceMxIn(any(localFaceMxIn,2),:);
tool = iFace1+length(localFaceMxIn)-1;
faceMx.in(iFace1:(tool),:) = localFaceMxIn;
faceMxSize.in = tool;
faceMx.sur(nS1:(nS1+nWalls-1),:) = localFaceMxSur;
faceMxSize.sur = (nS1+nWalls-1);

function [faceMx,  faceMxSize] =innerFaceMxMaker (~,faceMx,nWalls,csDensity,~,d,~,pIndex,upP,m,downP,j,k1,nBifBif,  faceMxSize, ptCoordMxSize)
n = faceMxSize.in;

nS1 = faceMxSize.sur;

iFace1 = n+1;
nS1 = nS1+1;
localFaceMxIn = zeros (d*2,6);
localFaceMxSur = zeros(nWalls,6);
iFace = 1;
nS = 1;
startFace = iFace;
startFace2 = nS;
offset = (d*(nBifBif-1));
upP(end+1,1:csDensity)=upP(1,1:csDensity);
CSupP(1,:) =[upP(1:(nWalls/4),csDensity+1)',upP(((nWalls)/4)+1,csDensity+1:end),flipud(upP(1:(nWalls/4),end))',fliplr(upP(1,csDensity+1:end-1))];
CSupP(2,:) =[upP(1:(nWalls/4),csDensity+1)',upP(((nWalls)/4),csDensity+1:end-1),flipud(upP(1:(nWalls/4),end-1))',fliplr(upP(1,csDensity+1:end-1)),0];
%%%%%%%%%%%%faceMx developmet%%%%%%%%%%%%%%%%%%%
if m == 1 && j == 1
    downP = pIndex;
    downP(end+1,1:csDensity) = downP(1,1:csDensity);
    CSdownP(1,:) = [pIndex(1:(nWalls/4),csDensity+1)',pIndex(((nWalls)/4)+1,csDensity+1:end),flipud(pIndex(1:(nWalls/4),end))',fliplr(pIndex(1,csDensity+1:end-1))];
    CSdownP(2,:) = [pIndex(1: (nWalls/4),csDensity+1)',pIndex(((nWalls)/4),csDensity+1:end-1),flipud(pIndex(1:(nWalls/4),end-1))',fliplr(pIndex(1,csDensity+1:end-1)),0];
else
    downP(end+1,1:csDensity) = downP(1,1:csDensity);
    CSdownP(1,:) = [downP(1:(nWalls/4),csDensity+1)',downP(((nWalls)/4)+1,csDensity+1:end),flipud(downP(1:(nWalls/4),end))',fliplr(downP(1,csDensity+1:end-1))];
    CSdownP(2,:) = [downP(1: (nWalls/4),csDensity+1)',downP(((nWalls)/4),csDensity+1:end-1),flipud(downP(1:(nWalls/4),end-1))',fliplr(downP(1,csDensity+1:end-1)),0];
end
for t = 1:csDensity
    i = 0;
    for a = 1:(nWalls)
        %%%alalin wall
        if t ==1
            i = i+1;
            localFaceMxIn(iFace,:) = [upP(i,t) upP(i+1,t) upP(i+1,t+1) upP(i,t+1) upP(i,t)+offset upP(i,t)+d+offset];
            localFaceMxSur(nS,:) = [upP(i,t) upP(i+1,t) downP(i+1,t) downP(i,t) 0 upP(i,t)+offset ];
            iFace = iFace + 1;
            nS = nS+1;
            if a ==1
                localFaceMxIn(iFace,:) = [upP(i,t) upP(i,t+1) downP(i,t+1) downP(i,t) upP(i,t)+offset upP(i,t)+nWalls-1+offset];
            else
                localFaceMxIn(iFace,:) = [upP(i,t) upP(i,t+1) downP(i,t+1) downP(i,t) upP(i,t)+offset upP(i,t)-1+offset];
            end
            iFace = iFace + 1;
        elseif t == csDensity
            i=i+1;
            localFaceMxIn(iFace,:) =  [upP(i,t) upP(i+1,t) CSupP(1,i+1) CSupP(1,i) upP(i,t)+offset upP(i,t)+d+offset];
            if a == 1
                localFaceMxIn(iFace+1,:) =  [upP(i,t) CSupP(1,i) CSdownP(1,i) downP(i,t) upP(i,t)+offset upP(i,t)+nWalls-1+offset];
            else
                localFaceMxIn(iFace+1,:) =  [upP(i,t) CSupP(1,i) CSdownP(1,i) downP(i,t) upP(i,t)+offset upP(i-1,t)+offset];
            end
            localFaceMxIn(iFace+2,:) =  [upP(i,t) upP(i+1,t) downP(i+1,t) downP(i,t) upP(i,t-1)+offset upP(i,t)+offset];
            localFaceMxIn(iFace+3,:) =  [CSupP(1,i) CSupP(1,i+1) CSdownP(1,i+1) CSdownP(1,i) upP(i,t)+offset CSupP(2,i)+offset];
            iFace = iFace + 4;
        elseif (t > 1 && t < csDensity)
            i = i+1;
            localFaceMxIn(iFace,:) = [upP(i,t) upP(i+1,t) upP(i+1,t+1) upP(i,t+1) upP(i,t)+offset upP(i,t)+d+offset];
            localFaceMxIn(iFace+1,:) = [upP(i,t) upP(i+1,t) downP(i+1,t) downP(i,t) upP(i,t-1)+offset upP(i,t)+offset ];
            if a ==1
                localFaceMxIn(iFace+2,:) = [upP(i,t) upP(i,t+1) downP(i,t+1) downP(i,t) upP(i,t)+offset upP(i,t)+nWalls-1+offset];
            else
                localFaceMxIn(iFace+2,:) = [upP(i,t) upP(i,t+1) downP(i,t+1) downP(i,t) upP(i,t)+offset upP(i,t)-1+offset];
            end
            iFace = iFace + 3;
        end
    end
end
for t = csDensity+1:csDensity+(nWalls/4)
    i=0;
    if  t ~= csDensity+(nWalls/4)
        for a = 1:(nWalls/4)
            if a == (nWalls/4)
                i=i+1;
                localFaceMxIn(iFace,:) = [upP(i,t) upP(i+1,t) upP(i+1,t+1) upP(i,t+1) upP(i,t)+offset upP(i,t)+d+offset];
                localFaceMxIn(iFace+1,:) = [upP(i,t+1) upP(i+1,t+1) downP(i+1,t+1) downP(i,t+1) upP(i,t)+offset upP(i,t+1)+offset];
                iFace = iFace+2;
            else
                i=i+1;
                localFaceMxIn(iFace,:) = [upP(i,t) upP(i+1,t) upP(i+1,t+1) upP(i,t+1) upP(i,t)+offset upP(i,t)+d+offset];
                localFaceMxIn(iFace+1,:) = [upP(i+1,t) upP(i+1,t+1) downP(i+1,t+1) downP(i+1,t) upP(i+1,t)+offset upP(i,t)+offset];
                localFaceMxIn(iFace+2,:) = [upP(i,t+1) upP(i+1,t+1) downP(i+1,t+1) downP(i,t+1) upP(i,t)+offset upP(i,t+1)+offset];
                iFace = iFace+3;
            end
        end
    else
        for a = 1:(nWalls/4)
            if a == (nWalls/4)
                i=i+1;
                localFaceMxIn(iFace,:) = [upP(i,t) upP(i+1,t) upP(i+1,t+1) upP(i,t+1) upP(i,t)+offset upP(i,t)+d+offset];
                iFace = iFace+1;
            else
                i=i+1;
                localFaceMxIn(iFace,:) = [upP(i,t) upP(i+1,t) upP(i+1,t+1) upP(i,t+1) upP(i,t)+offset upP(i,t)+d+offset];
                localFaceMxIn(iFace+1,:) = [upP(i+1,t) upP(i+1,t+1) downP(i+1,t+1) downP(i+1,t) upP(i+1,t)+offset upP(i,t)+offset];
                iFace = iFace+2;
            end
        end
    end
end
endFace = iFace-1;
endFace2 = nS-1;
if k1 == -1
    localFaceMxIn(startFace:endFace,1:6) = [localFaceMxIn(startFace:endFace,1:4),localFaceMxIn(startFace:endFace,6),localFaceMxIn(startFace:endFace,5)];
    localFaceMxSur(startFace2:endFace2,1:6) = [localFaceMxSur(startFace2:endFace2,1:4),localFaceMxSur(startFace2:endFace2,6),localFaceMxSur(startFace2:endFace2,5)];
end
localFaceMxIn = localFaceMxIn(any(localFaceMxIn,2),:);
tool = iFace1+length(localFaceMxIn)-1;
faceMx.in(iFace1:(tool),:) = localFaceMxIn;
faceMxSize.in = tool;
faceMx.sur(nS1:(nS1+nWalls-1),:) = localFaceMxSur;
faceMxSize.sur = (nS1+nWalls-1);

function [inMesh,squareP,ptCoordMx,indexMx,ptCoordMxSize] = computeInnermesh(points,center,wallDensity,fluidDensity,wall,fluid,nWalls,ptCoordMx,ptCoordMxSize)
d = 0:7;
totalP = (nWalls*(wallDensity + fluidDensity))+((nWalls/4)+1)^2;
tempPtCoord =[];
% tempPtCoord = zeros(totalP,3);
% tool = find(any(ptCoordMx,2),1,'last');
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
% toc
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

% toc

fluidP = cell(1,8);
for h=1:8
    fluidP{h} = zeros(nWalls/8,3,(fluidDensity+1));
    for u = 1:nWalls/8
        for c = 1:(fluidDensity+1)
            temp = dirSC{h}(u,:)*((c-1)/fluidDensity);
            fluidP{h}(u,:,c) = square{h}(u,:) + temp;
            %                  scatter3(fluidP{h}(u,1,c),fluidP{h}(u,2,c),fluidP{h}(u,3,c))
        end
        %                  squareDir2 = flipud(square{3})-square{1};
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
    %     tempPtCoord(1+((f-1)*nWalls):nWalls+((f-1)*nWalls),:) = inMesh{f};
    
    %                          for h =1:size(inMesh{f})
    %                         scatter3(inMesh{f}(h,1),inMesh{f}(h,2),inMesh{f}(h,3),'*r')
    %                          end
    inMesh{f} = [ptIndex,inMesh{f}(:,1:3)];
    indexMx =  [indexMx,ptIndex] ;
    ptCoordMxSize = ptCoordMxSize + nWalls;
end
for f = 1:(nWalls/4)+1
    ptIndex = ((ptCoordMxSize+1):(ptCoordMxSize+(nWalls/4)+1))';
    tempPtCoord = [tempPtCoord; squareP{f}];
    %         for h=1:size(squareP{f})
    %         scatter3(squareP{f}(h,1),squareP{f}(h,2),squareP{f}(h,3),'*g')
    %         end
    squareP{f} = [ptIndex,squareP{f}(:,1:3)];
    indexMx(1:(nWalls/4)+1,wallDensity+fluidDensity+f) = ptIndex ;
    ptCoordMxSize = ptCoordMxSize + (nWalls/4)+1;
end
ptCoordMx(ptCoordMxSize-totalP+1:(ptCoordMxSize),:) = tempPtCoord;
function startPoint2=computeStartPoint(point1,center1,velocity1,center2,velocity2,radius)
vector = point1-center1;
%tangent1 = normalize(c0_1-center1);
% scatter3(vector(:,1),vector(:,2),vector(:,3),'*k')
tangent1 = normalize(velocity1);
tangent2 = normalize(velocity2);
test = cross(vector,tangent1);
% if test ~= 0
normalVector1 = normalize(cross(vector,tangent1));
normalVector2 = normalize(cross(tangent2,normalVector1));
startPoint2 = center2+(normalVector2*radius);
% else
% startPoint2=ComputeProjectedPoint(point1,center1,center2,tangent2,radius);
% end
function projectedPoint=ComputeProjectedPoint(point1,centerPoint1,centerPoint2,tangent,radius1)
direction1=normalize(point1-centerPoint1);
projectedVector=direction1-((dot(direction1,normalize(tangent)))*normalize(tangent));
projectedPoint=(normalize(projectedVector)*radius1)+centerPoint2;
function Angle = computeAngle(point_1,middle_point,point_2)
axis1=point_1-middle_point;
axis2=point_2-middle_point;
Angle=atan2(norm(cross(axis1,axis2)),dot(axis1,axis2));

function velocity = computeVelocity(point1, c0_1, point2, c1_1, nSteps)
velocity = zeros(nSteps,3);
c0 = c0_1;
c1 = c1_1;
t = 0:(1/(nSteps)):1;
for i = (1:nSteps+1)
    T = t(i);
    velocity(i,:)=((3*(1-T)^2*(c0-point1))+(6*(T)*(1-T)*(c1-c0))+(3*(T)^2*(point2-c1)))';
end
function [Length] = computeArcLength(p0,p1,c0,c1, ~)
n=10;

t=0:(1/n):1;
bezier=zeros((n+1),3);
for i=1:(n+1)
    T=t(i);
    bezier(i,:)=((1-T)^3*p0+3*(T)*(1-T)^2*c0+3*(1-T)*(T)^2*c1+T^3*p1)';
end
for j=1:n
    Dist(j) = norm(bezier(j+1,:)-bezier(j,:))';
end
Length=sum(Dist);
function nSteps= computenSteps(arcLen,criteria,radius)
nSteps = round((arcLen/radius)*criteria);

if nSteps == 0 || nSteps == 1
    nSteps = 2;
end
function vector = normalize(aVector)
vector = aVector/norm(aVector);
function bezierCurve = computeBezier(point1, c0, point2, c1,nSteps,arcLen)
bezierCurve = zeros((nSteps+1),3);
t = 0:(1/(nSteps)):1;
for i = 1:(nSteps+1)
    T = t(i);
    bezierCurve(i,:)=((1-T)^3*point1+3*(T)*(1-T)^2*c0+3*(1-T)*(T)^2*c1+T^3*point2)';  
end
function distance=computeDistance(point1,point2)
x = [point1;point2];
distance = pdist(x);
















