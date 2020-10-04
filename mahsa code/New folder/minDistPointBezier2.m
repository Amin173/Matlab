 function [minPoint] = minDistPointBezier2(p0,c0,c1,p1,point)
% bezierDraw(p0,c0,c1,p1);
allP = [p0;c0;c1;p1];
% scatter3(allP(:,1),allP(:,2),allP(:,3),'filled', 'k');
% scatter3(point(:,1),point(:,2),point(:,3),'filled','m');
curve1 = allP;
for k = 1:5
    for i = 1:4
        distMatrix(i) = norm(point-curve1(i,:));
    end
    [~,col] = find(distMatrix == min(distMatrix(:)));
    if distMatrix(col) == 0
        minPoint = point;
    else
        if col == 1
            col = 2;
        elseif col == 4
            col = 3;
        end
        for i = 1:3
            lengthCurve1(i,:) = norm(curve1(i,:)-curve1(i+1,:));
        end
        sumLengthCurve1 = sum(lengthCurve1);
        if length(col) >= 2
            col = col(2);
            if col == 4
                col = 3;
            end
        end
        t1 = sum(lengthCurve1(1:(col-1)))/sumLengthCurve1;
        segment = findCloseSegment(p0,c0,c1,p1,t1,col,point);
        
        
        [curve1p0,curve1c0,curve1c1,curve1p1,splitPoint] = BezierSubdevision(curve1(1,:),curve1(2,:),curve1(3,:),curve1(4,:),t1);
        p0 = curve1p0(segment ,:);
        c0 = curve1c0(segment ,:);
        c1 = curve1c1(segment ,:);
        p1 = curve1p1(segment ,:);
        curveNext = [curve1p0(segment ,:);curve1c0(segment ,:);curve1c1(segment ,:);curve1p1(segment ,:)];
        curve1 = curveNext;
        %     scatter3(curve1(:,1),curve1(:,2),curve1(:,3))
        %     scatter3
        
        minPoint = splitPoint;
%         scatter3(minPoint(:,1),minPoint(:,2),minPoint(:,3),'filled','g')
    end
end
% scatter3(minPoint(:,1),minPoint(:,2),minPoint(:,3),'filled','g')

function [segment] = findCloseSegment(p0,c0,c1,p1,t,col,point)
if col == 2 || col == 3
    t_1 = t-0.05;
    t_2 = t+0.05;
    testPoint1 = ((1-t_1)^3.*p0+3*(t_1)*(1-t_1)^2.*c0+3*(1-t_1)*(t_1)^2.*c1+t_1^3.*p1);
    testPoint2 = ((1-t_2)^3.*p0+3*(t_2)*(1-t_2)^2.*c0+3*(1-t_2)*(t_2)^2.*c1+t_2^3.*p1);
    %     scatter3(testPoint1(:,1),testPoint1(:,2),testPoint1(:,3),'k')
    %     scatter3(testPoint2(:,1),testPoint2(:,2),testPoint2(:,3),'b')
    norm1 = pdist([testPoint1;point]);
    norm2 = pdist([testPoint2;point]);
    if norm1 < norm2
        %         col_2 = col+1;
        segment = 1;
    else
        %         col_2 = col-1;
        %         1;
        segment = 2;
    end
else
end
%     if col == 1
%     col_2 = col+1;
%     else
%     col_2 = 3;
%     end
% end
% sortCol = sort([col , col_2]);
% segment = min([col , col_2]);
