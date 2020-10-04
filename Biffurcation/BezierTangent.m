function velocityVec = BezierTangent(t,splineCtrlPts)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
b0 = 3*(1-t).^2;
b1 = 6*(1-t).*t;
b2 = 3*t.^2;

velocityVec = b0*(splineCtrlPts(:,2)-splineCtrlPts(:,1))+...
    b1*(splineCtrlPts(:,3)-splineCtrlPts(:,2))+...
    b2*(splineCtrlPts(:,4)-splineCtrlPts(:,3));
end

