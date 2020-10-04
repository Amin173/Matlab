function d=bmin(obs,p1,p2)
obs=[obs obs(:,1)];
for i=1:4
    v1=obs(:,i);
    v2=obs(:,i+1);
      [B,V]= DistBetween2Segment(p1,p2,v1,v2)
      D(i)=B;
      C(:,i)=V;
end
d=min(D);
c=C(:,D==d);
d=[d;0];
c=c(:,1);
d=[d c];
end