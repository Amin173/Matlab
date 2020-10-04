clc
clear
close all
%% 
% Generating data points
x1=(0:0.01:1)';
x2=(0:0.01:1)';
B1=sin((x1)*pi);
B2=-sin((x2)*pi);
plot(x1,B1,x2,B2);
% Bezier coefficients
t=setTValueForPts([B1,x],length(x)); %the same for both splines
b0 = (1-t).^3;
b1 = 3*t.*(1-t).^2;
b2 = 3*(1-t).*t.^2;
b3 = t.^3;
% Obtaining control points
A=[x;B1;x;flipud(B2)];
b=[b3/2 b0/2 b0/2+b1 b3/2+b2 zeros(length(t),4);zeros(length(t),4) b3/2 b0/2 b0/2+b1 b3/2+b2;b0/2+b1 b2+b3/2 b3/2 b0/2 zeros(length(t),4); zeros(length(t),4) b0/2+b1 b2+b3/2 b3/2 b0/2];
c=b\A;
% Reconstructing data points by bezier function
c=reshape(c,4,2);
c=c';
B_rep_1=b0'.*(c(:,1)+c(:,4))/2+b1'.*c(:,1)+b2'.*c(:,2)+b3'.*(c(:,2)+c(:,3))/2;
B_rep_2=b0'.*(c(:,1)+c(:,4))/2+b1'.*c(:,3)+b2'.*c(:,4)+b3'.*(c(:,2)+c(:,3))/2;
hold on
plot(B_rep_1(1,:),B_rep_1(2,:),':')
plot(B_rep_2(1,:),B_rep_2(2,:),':')
scatter(c(1,1),c(2,1),'k');  scatter(c(1,2),c(2,2),'b');  scatter(c(1,3),c(2,3),'m');  scatter(c(1,4),c(2,4),'g');
text(c(1,1),c(2,1),'C1')
text(c(1,2),c(2,2),'C2')
text(c(1,3),c(2,3),'C3')
text(c(1,4),c(2,4),'C4')
title('My solution')

%%  Grant solution
% find the bezier curves
SolutionLoop = 1;
% for SolutionLoop = 1:nSplines
%first bezier
xPos1 = x;
yPos1 = B1;
t = setTValueForPts([B1,x],length(x));
%     tMatrix = zeros(length(xPos1),4);
tMatHold = zeros(length(xPos1),4);
% set up the matrix for calculating bezier curves --> has to be re-written every time a bezier curve is added
for i = 1:length(xPos1)
    % p1, p2, p3, p4 (p0 and p5 are assumed later)
    tMatrixDualBezierTop(i,:)= [0.5*t(i)^3, 0.5*(1-t(i))^3, 0.5*(1-t(i))^3+3*t(i)*(1-t(i))^2, 0.5*t(i)^3+3*(1-t(i))*t(i)^2];
    tMatHold(i,:)= [(1-t(i))^3 3*t(i)*(1-t(i))^2, 3*(1-t(i))*t(i)^2 t(i)^3];   %for reconstruction later
end
% second bezier
xPos2 = x;
yPos2 = B2;
t2 = setTValueForPts([flipud(B2),x],length(x));
%     xPos2 = xPoints{2*SolutionLoop};
%     yPos2 = yPoints{2*SolutionLoop};
%     t2 = tTotal{2*SolutionLoop};
%     tMatrix2 = zeros(length(xPos1),4);
tMatHold2 = zeros(length(xPos1),4);
for i = 1:length(xPos2)
    %         tMatrixDualBezierBottom(i,:)= [0.5*(1-t2(i))^3+3*t2(i)*(1-t2(i))^2, 0.5*t2(i)^3+3*t2(i)^2*(1-t2(i)), 0.5*t2(i)^3, 0.5*(1-t2(i))^3];
    tMatrixDualBezierBottom(i,:)= [0.5*(1-t2(i))^3+3*t2(i)*(1-t2(i))^2, 0.5*t2(i)^3+3*t2(i)^2*(1-t2(i)), 0.5*t2(i)^3, 0.5*(1-t2(i))^3];
    tMatHold2(i,:)= [(1-t2(i))^3 3*t2(i)*(1-t2(i))^2, 3*(1-t2(i))*t2(i)^2 t2(i)^3];   %for reconstruction later
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
figure(2)
plot(xPos1,yPos1,xPos2,yPos2,tMatHold*PVect1X',tMatHold*PVect1Y',tMatHold2*PVect2X',tMatHold2*PVect2Y');
hold on
scatter(pX(1),pY(1),'k');  scatter(pX(2),pY(2),'b');  scatter(pX(3),pY(3),'m');  scatter(pX(4),pY(4),'g');
%     scatter(PVect1X,PVect1Y,'k');
%     scatter(PVect2X,PVect2Y,'m');
%     plot(xPos1,yPos1,xPos2,yPos2,tMatHold*PVect1X',tMatHold*PVect1Y',tMatHold2*PVect2X',tMatHold2*PVect2Y');
text(pX(1),pY(1),'C1')
text(pX(2),pY(2),'C2')
text(pX(3),pY(3),'C3')
text(pX(4),pY(4),'C4')
title('Your solution')


%%
function T = setTValueForPts(dataSplitIntoSplines,nPtsInSpline)
splinePts = dataSplitIntoSplines;
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
T = TLocal;

end