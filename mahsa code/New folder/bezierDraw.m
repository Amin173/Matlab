function drawBezier = bezierDraw(p0,c0,c1,p1,style)
for i = 1:size(p0,1)
                t = 0:0.05:1;
                bezierCurve = zeros(length(t),3);
                for w = 1:(length(t))
                    T = t(w);
                    bezierCurve(w,:)=((1-T)^3*p0(i,:)+3*(T)*(1-T)^2*c0(i,:)+3*(1-T)*(T)^2*c1(i,:)+T^3*p1(i,:));
                end
                hold on
                if style == 1
         
                    plot3(bezierCurve(:,1),bezierCurve(:,2),bezierCurve(:,3),'b.-','linewidth',1.5);
                else
                    plot3(bezierCurve(:,1),bezierCurve(:,2),bezierCurve(:,3),'color',[rand() rand() rand()],'linewidth',1.5);
                end
                
end
         axis equal