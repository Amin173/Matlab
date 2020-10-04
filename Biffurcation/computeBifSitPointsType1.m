
function bifurSitPoints=computeBifSitPointsType1(point1,point2,bifCoord,wallDistribution)
counter= 0;
temp1 = point1-bifCoord;
tempAbs1 = (temp1(1)^2+temp1(2)^2+temp1(3)^2)^.5;

temp2 = point2-bifCoord;
tempAbs2 = (temp2(1)^2+temp2(2)^2+temp2(3)^2)^.5;
r = zeros (wallDistribution,1);
if tempAbs1 >= tempAbs2
    for t= 1:-(1/(wallDistribution)):0
        counter= counter+1;
        r(counter,:)= [t 1-t]*[tempAbs1;tempAbs2];
    end
else
    for t= 1:-(1/(wallDistribution)):0
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
for h= 1:(wallDistribution)
    teta= (h-1)*(temp/wallDistribution);
    counter1= counter1+1;
    bifurSitPoints(h,1:3)=bifCoord + (cos(teta)*r(counter1)*b) + (r(counter1)*sin(teta)*b2);
    
end
end