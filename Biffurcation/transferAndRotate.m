function newCrossSection = transferAndRotate(crossSection,J,ppi1,ppi2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
crossSection = reshape(crossSection,3,[]);
newCrossSection = zeros(size(crossSection));
for i = 1:size(crossSection,2)
    newCrossSection(:,i) = J*(crossSection(:,i)-ppi1)+ppi2;
end
end

