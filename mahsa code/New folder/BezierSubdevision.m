%%%bezier subdevision algorithm 2.19.2015 Mahsa Ghaffari
function [p0,c0,c1,p1,finalPoint] = BezierSubdevision(p0,c0,c1,p1,t)
allP = [p0;c0;c1;p1];
% scatter3(allP(:,1),allP(:,2),allP(:,3),'filled', 'k')
% bezierDraw(p0,c0,c1,p1);
Dir = zeros(3,3);
Dir(1,:) = c0 - p0;
Dir(2,:) = c1 - c0;
Dir(3,:) = p1 - c1;
pointsDir = t.*Dir;
points = [p0;c0;c1]+pointsDir;
% scatter3(points(:,1),points(:,2),points(:,3),'*k')

Dir2 = t*(points(2:3,:) - points(1:2,:));
points2 = points(1:2,:) + Dir2;
% scatter3(points2(:,1),points2(:,2),points2(:,3),'*b');

finalPoint = t*(points2(2,:)-points2(1,:))+ points2(1,:);
%  scatter3(finalPoint(:,1),finalPoint(:,2),finalPoint(:,3),'filled');
p0 = [p0; finalPoint];
c0 = [points(1,:); points2(2,:)];
c1 = [points2(1,:);points(end,:)];
p1 = [finalPoint;p1];

  

