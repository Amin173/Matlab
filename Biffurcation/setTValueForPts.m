function T = setTValueForPts(dataSplitIntoSplines,nPtsInSpline)
%This function computes the t values for a given curve in three
%dimentional space

splinePts = dataSplitIntoSplines;
totalSubSplineLength = 0;
TLocal = zeros(nPtsInSpline,1);
%     get total length of spline
for iPtInSpline = 1:nPtsInSpline-1
    p1 = splinePts(iPtInSpline,:);
    p2 = splinePts(iPtInSpline+1,:);
    dist(iPtInSpline) = sqrt((p2(1)-p1(1))^2 + (p2(2)-p1(2))^2+ (p2(3)-p1(3))^2);
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