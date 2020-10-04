function HW3_doubleBezier_cleaned
tt=[-0.5:0.01:0.5];
data1=-sin(tt*pi);
data2=sin(tt*pi);
% data1 = load('Pic 1 Left.csv');
% data2 = load('Pic 1 Right.csv');
% data3 = load('Pic 2 Left.csv');
% data4 = load('Pic 2 Right.csv');
% data4 = flipud(data4);
% dataTemp = {[data1; flipud(data2)], [data3; flipud(data4)]};
imageData1 = [data1; flipud(data2)];
% imageData2 = [data3; flipud(data4)];
bezierPieces = 2;
% imageNum = 2;
nSplinePoints = length(imageData1)/bezierPieces;
% size(2) = length(dataTemp{2})/bezierPieces;
processImage(imageData1,nSplinePoints,2);
end

function processImage(imageData,nSplinePoints,nSplines)
%% split the data for each bezier piece
dataSplitIntoPieces = splitDataIntoNPieces(imageData,nSplines,nSplinePoints);
T = setTValueForPts(dataSplitIntoPieces,nSplinePoints,nSplines);
%% find the bezier curves
SolutionLoop = 1;
% for SolutionLoop = 1:nSplines
%first bezier
data = dataSplitIntoPieces{SolutionLoop};
xPos1 = data(:,1);
yPos1 = data(:,2);
t = T{SolutionLoop};
%     tMatrix = zeros(length(xPos1),4);
tMatHold = zeros(length(xPos1),4);
%% set up the matrix for calculating bezier curves --> has to be re-written every time a bezier curve is added
for i = 1:length(xPos1)
    % p1, p2, p3, p4 (p0 and p5 are assumed later)
    tMatrixDualBezierTop(i,:)= [0.5*t(i)^3, 0.5*(1-t(i))^3, 0.5*(1-t(i))^3+3*t(i)*(1-t(i))^2, 0.5*t(i)^3+3*(1-t(i))*t(i)^2];
    tMatHold(i,:)= [(1-t(i))^3, 3*t(i)*(1-t(i))^2, 3*(1-t(i))*t(i)^2, t(i)^3];   %for reconstruction later
end
% second bezier
data = dataSplitIntoPieces{2};
xPos2 = data(:,1);
yPos2 = data(:,2);
t2 = T{2};
%     xPos2 = xPoints{2*SolutionLoop};
%     yPos2 = yPoints{2*SolutionLoop};
%     t2 = tTotal{2*SolutionLoop};
%     tMatrix2 = zeros(length(xPos1),4);
tMatHold2 = zeros(length(xPos1),4);
for i = 1:length(xPos2)
    %         tMatrixDualBezierBottom(i,:)= [0.5*(1-t2(i))^3+3*t2(i)*(1-t2(i))^2, 0.5*t2(i)^3+3*t2(i)^2*(1-t2(i)), 0.5*t2(i)^3, 0.5*(1-t2(i))^3];
    tMatrixDualBezierBottom(i,:)= [0.5*(1-t2(i))^3+3*t2(i)*(1-t2(i))^2, 0.5*t2(i)^3+3*t2(i)^2*(1-t2(i)), 0.5*t2(i)^3, 0.5*(1-t2(i))^3];
    tMatHold2(i,:)= [(1-t2(i))^3, 3*t2(i)*(1-t2(i))^2, 3*(1-t2(i))*t2(i)^2, t2(i)^3];   %for reconstruction later
end
tMatDualBezier = [tMatrixDualBezierTop; tMatrixDualBezierBottom];
bVectSingleBezierX = [xPos1; xPos2];
bVectSingleBezierY = [yPos1; yPos2];

%     tNewTotal{SolutionLoop} = tMatrix;
pX = tMatDualBezier\bVectSingleBezierX;
pY = tMatDualBezier\bVectSingleBezierY;
PVect1X = [(pX(1)+pX(4))/2 pX(1) pX(2) (pX(2)+pX(3))/2];
PVect1Y = [(pY(1)+pY(4))/2 pY(1) pY(2) (pY(2)+pY(3))/2];
PVect2X = [(pX(2)+pX(3))/2 pX(3) pX(4) (pX(1)+pX(4))/2];
PVect2Y = [(pY(2)+pY(3))/2 pY(3) pY(4) (pY(1)+pY(4))/2];

% end
plot(xPos1,yPos1,xPos2,yPos2,tMatHold*PVect1X',tMatHold*PVect1Y',tMatHold2*PVect2X',tMatHold2*PVect2Y');
hold on
scatter(pX(1),pY(1),'k');  scatter(pX(2),pY(2),'b');  scatter(pX(3),pY(3),'m');  scatter(pX(4),pY(4),'g');
%     scatter(PVect1X,PVect1Y,'k');
%     scatter(PVect2X,PVect2Y,'m');
%     plot(xPos1,yPos1,xPos2,yPos2,tMatHold*PVect1X',tMatHold*PVect1Y',tMatHold2*PVect2X',tMatHold2*PVect2Y');
end

function T = setTValueForPts(dataSplitIntoSplines,nPtsInSpline,nSplines)
for iSpline = 1:nSplines
    splinePts = dataSplitIntoSplines{iSpline};
    totalSubSplineLength = 0;
    TLocal = zeros(nPtsInSpline,1);
    %     get total length of spline
    for iPtInSpline = 1:nPtsInSpline-1
        p1 = splinePts(iPtInSpline,:);
        p2 = splinePts(iPtInSpline+1,:);
        dist(iPtInSpline) = sqrt((p2(1)-p1(1))^2 + (p2(2)-p1(2))^2);
        totalSubSplineLength = totalSubSplineLength + dist(iPtInSpline);
    end
    %     get t value of each point
    accumDist = 0;
    TLocal(1) = 0;
    for iPtInSpline = 1:nPtsInSpline-1
        distLocal = dist(iPtInSpline);
        accumDist = accumDist + distLocal;
        TLocal(iPtInSpline+1) = accumDist/totalSubSplineLength;
    end
    T{iSpline} = TLocal;
end
end

function dataSplitIntoPieces = splitDataIntoNPieces(imageData,nSplines,nSplinePoints)
for i = 1:nSplines
    dataSubset = imageData(1+(i-1)*nSplinePoints:i*nSplinePoints,:);
    dataSplitIntoPieces{i} = dataSubset;
end
end