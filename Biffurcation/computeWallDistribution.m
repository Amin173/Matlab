function wallDistribution = computeWallDistribution(controlPoints,separationPoints,bifCoord,nWalls)
%%the wallDistribution1 shows distribution of walls around bifucation
%%site,the distribution start from bifurcation axis in clock wise each
%%column of wallDistribution shows the wallDistribution from C1 in clock
%%wise for each vessels
halfWalls=nWalls/2;
criterion=pi/halfWalls;

wallDistribution(1)=round((computeAngle(separationPoints(1,:),bifCoord,controlPoints(1,:)))/criterion);
wallDistribution(2)=round((computeAngle(separationPoints(1,:),bifCoord,controlPoints(2,:)))/criterion);
wallDistribution(3)=round((computeAngle(separationPoints(2,:),bifCoord,controlPoints(2,:)))/criterion);
wallDistribution(4)=round((computeAngle(separationPoints(2,:),bifCoord,controlPoints(1,:)))/criterion);
wallDistribution(5)=round((computeAngle(separationPoints(3,:),bifCoord,controlPoints(2,:)))/criterion);
wallDistribution(6)=round((computeAngle(separationPoints(3,:),bifCoord,controlPoints(1,:)))/criterion);

if wallDistribution(1) == 0
    wallDistribution(1)= 1;
    wallDistribution(2)= wallDistribution(2)-1;
end

if wallDistribution(3) == 0
    wallDistribution(3)= 1;
    wallDistribution(4)= wallDistribution(4)-1;
end

if wallDistribution(6) == 0
    wallDistribution(6)= 1;
    wallDistribution(5)= wallDistribution(5)-1;
end

wallDistribution1=[wallDistribution(1);wallDistribution(2);wallDistribution(3);wallDistribution(4)];
wallDistribution2=[wallDistribution(1);wallDistribution(2);wallDistribution(5);wallDistribution(6)];
wallDistribution=cat(2,wallDistribution1,wallDistribution2,wallDistribution1);

end