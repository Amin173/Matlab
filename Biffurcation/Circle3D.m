function points=Circle3D(radius,center,normal)
theta=0:pi/20:2*pi;
center=reshape(center,1,3);
normal=reshape(normal,1,3);
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
hold on
plot3(points(1,:),points(2,:),points(3,:),'b-','LineWidth',2);

end