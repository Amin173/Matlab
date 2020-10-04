function transGeom = SurfMesh(t1,t2,NSegments,e1,e2,e3,crossSection,splineCtrlPts)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dt=(t2-t1)/NSegments;
e1 = reshape(e1,[3,1]);
e2 = reshape(e2,[3,1]);
e3 = reshape(e3,[3,1]);
e11=e1;
e12=e2;
e13=e3;
tcutoff=0.5;
saveCrossSection{1}=crossSection;
ez{1}=e3;
tt=t1:dt:t2;
ppi(:,1)=splineCtrlPts(:,1);
A=crossSection;
for i=2:NSegments
    ppi(:,i)=BezierSpline(tt(i),splineCtrlPts);
    vi=BezierTangent(tt(i),splineCtrlPts);
%         if tt<tcutoff
%             alfa = tt/tcutoff;
%         elseif tt>(1-tcutoff)
%             alfa = (1-tt)/tcutoff;
%         else
%             alfa = 1;
%         end
    alfa=(tt(i)/tcutoff)*(tt(i)<tcutoff)+(1-tt(i))/tcutoff*(tt(i)>(1-tcutoff))+...
        1*(tt(i)>tcutoff&tt(i)<(1-tcutoff));
    nn1 = e13*(1-alfa)+alfa*vi;
    aaProj = dot(nn1,e11)*normalize(nn1);
    e1i = normalize(e11-aaProj);
    e3i = normalize(nn1);
    e2i = normalize(cross(e3i,e1i));
    JJ = [e1i,e2i,e3i];
    R=JJ/[e11,e12,e13];
    ez{i}=e3i;
    A = transferAndRotate(A,R,ppi(:,i-1),ppi(:,i));
    saveCrossSection{i} = A;
    e11 = e1i; 
    e12 = e2i; 
    e13 = e3i;


end
transGeom.SMesh=saveCrossSection;
transGeom.Nvector=ez;
end

