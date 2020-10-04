function newCrossSection = fixRoll(A,minRoll,V,ii,jj)
%normal vetor of cross section plane (e3)
center=mean(A,2);
angle=180*minRoll/(pi);
nv=V/norm(V);
angle=angle*(ii>jj)-angle*(jj>=ii);
for i=1:size(A,2)
    r=A(:,i)-center;
    r2=r/norm(r);
r2=rotVecAroundArbAxis(r2',nv',angle);
newCrossSection(:,i)=norm(r)*r2'+center;
end
end

