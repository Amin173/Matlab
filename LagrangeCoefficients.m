function [l] = LagrangeCoefficients(x,points)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% points=reshape(points,[],2);
% x=reshape(x,[1,2]);
n=length(points);
l=ones(n,1);
for i=1:n
    for j=1:n
        if i~=j
    l(i)=(x-points(j))*l(i);
    l(i)=l(i)/(points(i)-points(j));
        end
    end
end

