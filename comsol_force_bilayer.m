Read=csvread('Untitled.csv',5,0);
g=Read(:,1);
S184=Read(:,2);
force=Read(:,3);
i=0;
for j=5:5:40
    i=i+1;
    F(:,i)=force(g==j);
    S(:,i)=S184(g==j);
end
for j=1:i
    hold on
plot(S(:,j),F(:,j),'LineWidth',2)
end
legend('e2=5','e2=10','e2=15','e2=20','e2=25','e2=30','e2=35','e2=40')
for j=1:i
    hold on
plot(S(F(:,j)==max(F(:,j)),j),max(F(:,j)),'k+')
end

grid
grid minor
set(gca,'Fontname','Times New Roman','Fontsize',24)
xlabel('Dielectric constant of the dielectric layer')
ylabel('Electrostatic Adhesion Pressure (kPa)')