function XY=xy(a,b)
global l1 lc2 lc1 l2
XY(:,1)=l1*cos(a)+l2*cos(a+b);
XY(:,2)=l1*sin(a)+l2*sin(a+b);
end