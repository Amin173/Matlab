function d = rho(pt,obs)

obs=[obs obs(:,1)];
for i=1:4

 x =pt'; %some point
a = obs(:,i)'; %segment points a,b
b = obs(:,i+1)';

d_ab = norm(a-b);
d_ax = norm(a-x);
d_bx = norm(b-x);

if dot(a-b,x-b)*dot(b-a,x-a)>=0
    A = [a,1;b,1;x,1];
    dist = abs(det(A))/d_ab;        
else
    dist = min(d_ax, d_bx);
end
D(i) = dist;
%this is equivalent to the following line for a single point
%distance=norm(cross(v1-v2,pt-v2))/norm(v1-v2)
end
d=min(D);
end
