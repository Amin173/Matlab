function [ex,ey,ez] = bifaxis(ctrlpt,separationpt,bifpt)
%bifaxis This function generates a right handed coordinate axis at the
%input given biffurcation point (bifpt), which encloses the given seperation point
%as input (seperationpt.  
%   The x axis connects the biffurcation point to the control point given
%   as input.
A = separationpt-bifpt;
eA=A/norm(A);
ex = ctrlpt-bifpt;
ex=ex/norm(ex);

ex = ex*(dot(eA,ex)>=0)-ex*(dot(eA,ex)<0);
ez=cross(eA,ex);

if(norm(ez)~=0)
ez = ez*(dot(eA,ez)>=0)-ez*(dot(eA,ez)<0);
ez=ez/norm(ez);
ey = cross(ez,ex);
ey=ey/norm(ey);
if (dot(eA,ey)<0)
    ey=ez;
    ez=cross(ex,ey);
    ez=ez/norm(ez);
end
end
ex=ex+bifpt;
ey=ey+bifpt;
ez=ez+bifpt;
end
