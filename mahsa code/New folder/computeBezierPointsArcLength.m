function points = computeBezierPointsArcLength(p0,c0,c1,p1,t)
% stepsize = 1/(n-1);
% t = 0.0;
n = length(t);
points = zeros(n,3);
% for i = 1:n    
    points(:,1) = p0(1)*((1-t).^3) + 3*c0(1).*t.*(1-t).^2 + 3*c1(1).*t.^2.*(1-t) + p1(1).*t.^3;
    points(:,2) = p0(2)*((1-t).^3) + 3*c0(2).*t.*(1-t).^2 + 3*c1(2).*t.^2.*(1-t) + p1(2).*t.^3;
    points(:,3) = p0(3)*((1-t).^3) + 3*c0(3).*t.*(1-t).^2 + 3*c1(3).*t.^2.*(1-t) + p1(3).*t.^3;
%     t = t+stepsize;
% end