% clc
% clear 
% close all
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
function [x,y,z] = genData(nSplines,nSplinePts)
t = 0: 1/nSplinePts:1;

b0 = (1-t).^3;
b1 = 3*t.*(1-t).^2;
b2 = 3*(1-t).*t.^2;
b3 = t.^3;
%% the second amount of control point which are chosen by using continuity
P0 = [0.1,0.2,0.1; 1.1,1.2,1.1; 2.1,2.1,2.1; 2.1,2.1,2.1; 1.1,1.2,1.1];
P1 = [0.3,0.4,1; 1.3,1.3,1.1; 2.3, 2.5,2.8; 2.3,2.3,2.8; 1.5,1.5,-1.8];
P2 = [0.6,0.7,1; 1.6,1.5,1.1; 2.6,3.3,3.5; 3.6,2.6,4; 1.7,1.9,-2.2];
P3 = [1.1,1.2,1.1; 2.1,2.1,2.1; 3,4,5; 4,3,5; 3,2,-2];
for i=1:size(P0,1)
x(:,i) = P0(i,1)*b0+P1(i,1)*b1+P2(i,1)*b2+P3(i,1)*b3;
y(:,i) = P0(i,2)*b0+P1(i,2)*b1+P2(i,2)*b2+P3(i,2)*b3;
z(:,i) = P0(i,3)*b0+P1(i,3)*b1+P2(i,3)*b2+P3(i,3)*b3;
end
x=reshape(x,nSplinePts+1,nSplines);
y=reshape(y,nSplinePts+1,nSplines);
z = reshape(z,nSplinePts+1,nSplines);

end
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
% 
% 
