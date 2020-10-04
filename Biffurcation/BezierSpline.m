function BezierPt = BezierSpline(t,splineCtrlPts)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
b0 = (1-t).^3;
b1 = 3*t.*(1-t).^2;
b2 = 3*(1-t).*t.^2;
b3 = t.^3; 
BezierPt = splineCtrlPts(:,1)*b0+splineCtrlPts(:,2)*b1+...
    splineCtrlPts(:,3)*b2+splineCtrlPts(:,4)*b3;
end

