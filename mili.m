clc
clear 
close all
%% plot the lines and endpoints of branch
% pts = [0,0; 1,1; 2,2 ; 4,3; 2,0.5; 4,2];
% faces = [pts(1,:),pts(2,:);pts(2,:),pts(3,:);pts(3,:),pts(4,:);pts(3,:),pts(6,:);pts(2,:),pts(5,:)];
% 
% for i = 1 : length(faces)
%  plot([faces(i,1),faces(i,3)], [faces(i,2),faces(i,4)],'k', 'LineWidth', 2);
%  hold on
% end

% for i = 1 : length(pts)
% plot(pts(i,1), pts(i,2),'ko');
% end
% hold off
%% plot the bezier curves of each line (both control points in one side of the curve)
t = 0: .001:1;

b0 = (1-t).^3;
b1 = 3*t.*(1-t).^2;
b2 = 3*(1-t).*t.^2;
b3 = t.^3;
%% the first amount of the control points which are chosen randomly
% P0 = [0,0; 1,1 ; 2,2; 2,2; 1,1];
% P1 = [0.3,0.4; 1.3,1.5; 2.6, 2.8; 2.5,1.8; 1.2,0.3];
% P2 = [0.6,0.8; 1.6,1.8; 3.6,2.8; 3.5,1.8; 1.8,0.4];
% P3 = [1,1; 2,2; 4,3; 4,2; 2,0.5];
%% the second amount of control point which are chosen by using continuity
% P0 = [0,0; 1,1 ; 2,2; 2,2; 1,1];
% P1 = [-0.5,0.4; 1.3,1.2; 2.6, 2.8; 2.5,1.8; 1.2,0.3];
% P2 = [0.7,1; 1.6,1.8; 3.6,2.8; 3.5,1.8; 1.8,0.4];
% P3 = [1,1; 2,2; 4,3; 4,2; 2,0.5];
% 
% 
% for i=1:5
% x = P0(i,1)*b0+P1(i,1)*b1+P2(i,1)*b2+P3(i,1)*b3;
% y = P0(i,2)*b0+P1(i,2)*b1+P2(i,2)*b2+P3(i,2)*b3;
% 
% figure(1)
% plot(x,y)
% hold on
% end
% hold off
%% plot the bezier curves of each line (control points in different sides of the curve)
% t = 0: .001:1;
% 
% b0 = (1-t).^3;
% b1 = 3*t.*(1-t).^2;
% b2 = 3*(1-t).*t.^2;
% b3 = t.^3;
% 
% P0 = [0,0; 1,1 ; 2,2; 2,2; 1,1];
% P1 = [0.2,0.8; 1.2,1.8; 2.8, 3.5; 2.5,1.2; 2.5,1.2];
% P2 = [0.8,0.2; 1.8,1.2; 3.5,2.2; 3.5,2.5; 1.5,1];
% P3 = [1,1; 2,2; 4,3; 4,2; 2,0.5];
% 
% for i=1:5
% x = P0(i,1)*b0+P1(i,1)*b1+P2(i,1)*b2+P3(i,1)*b3;
% y = P0(i,2)*b0+P1(i,2)*b1+P2(i,2)*b2+P3(i,2)*b3;
% 
% figure(1)
% plot(x,y)
% hold on
% end
% hold off

%% finding the best Bezier curve for each line
% t = 0: .001:1;
% 
% b0 = (1-t).^3;
% b1 = 3*t.*(1-t).^2;
% b2 = 3*(1-t).*t.^2;
% b3 = t.^3;
% 
% P0 = [0,0; 1,1 ; 2,2; 2,2; 1,1];
% P1 = P1(5,2);
% P2 = P2(5,2);
% P3 = [1,1; 2,2; 4,3; 4,2; 2,0.5];
% 
% for i = 1:4;
% % P1 = zeros(5,2);
% % P2 = zeros(5,2);
% 
% %  P0(i+1,1) = 0.5.* (P2(i,1)+P1(i+1,1))   
% %  P0(i+1,2) =0.5.*(P2(i,2)+P1(i+1,2))   
% 
% 2.*P0(i+1,1)==P1(i+1,1)-P2(i,1);
% -2.*P1(i+1,1)+P2(i+1,1)-P1(i,1)+2.*P2(i,1)==0;
% 
% 2.*P0(i+1,2)-P1(i+1,2)-P2(i,2)==0;
% -2.*P1(i+1,2)+P2(i+1,2)-P1(i,2)+2.*P2(i,2)==0;
% 
% P0(i,1)-2.*P1(i,1)+P2(i,1)==0;
% P0(i,2)-2.*P1(i,2)+P2(i,2)==0;
% 
% P1(i,1)-2.*P2(i,1)+P0(i+1,1)==0;
% P1(i,2)-2.*P2(i,2)+P0(i+1,2)==0;
% 
% end
% 
% 
% for i=1:5
% x = P0(i,1)*b0+P1(i,1)*b1+P2(i,1)*b2+P3(i,1)*b3;
% y = P0(i,2)*b0+P1(i,2)*b1+P2(i,2)*b2+P3(i,2)*b3;
% 
% figure(1)
% plot(x,y)
% hold on
% end
% hold off
%% 3D
P0 = [0,0,0; 1,1,1 ; 2,2,2; 2,2,2; 1,1,1];
P1 = [0.3,0.3,0.7; 1.3,1.3,1.5; 2.3, 2.5,2.8; 2.3,2.3,2.8; 1.3,1.3,2.5];
P2 = [0.6,0.6,1.2; 1.6,1.8,1.8; 2.6,3.3,3.5; 3.6,2.6,4.8; 1.6,1.6,4.2];
P3 = [1,1,1; 2,2,2; 3,4,5; 4,3,5; 3,1,7];


for i=1:5
x = P0(i,1)*b0+P1(i,1)*b1+P2(i,1)*b2+P3(i,1)*b3;
y = P0(i,2)*b0+P1(i,2)*b1+P2(i,2)*b2+P3(i,2)*b3;
z = P0(i,3)*b0+P1(i,2)*b1+P2(i,2)*b2+P3(i,2)*b3;

figure(1)
plot3(x,y,z)
hold on
end
hold off

Rc=0.3;
r = Rc;
[x,y,z] = tubeplot(x,y,z,r,10,100);
surf(X,Y,Z,'EdgeColor','none','FaceAlpha',0.01,'FaceColor','b')