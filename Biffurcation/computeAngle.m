function Angle = computeAngle(point_1,middle_point,point_2)
axis1=point_1-middle_point;
axis2=point_2-middle_point;
Angle=atan2(norm(cross(axis1,axis2)),dot(axis1,axis2));
end