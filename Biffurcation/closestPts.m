function [r1,c,r2,ii,jj] = closestPts(firstCSection,secondCSection)
for i=1:size(firstCSection,2)
   r(:,i)=sum((secondCSection-firstCSection(:,i)).^2)';
end
ii=1:length(min(r)); 
ii=ii(min(r)==min(min(r)));
jj=1:length(r(:,ii));
jj=jj(r(:,ii)'==min(r(:,ii)));
%% defining the closest vectors
r1=firstCSection(:,ii);
r2=secondCSection(:,jj);
center1=mean(firstCSection,2);
center2=mean(secondCSection,2);
c = (center1+center2)/2;
end

