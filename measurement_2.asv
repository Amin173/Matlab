function [z mew] = measurement_2(mew0,U,Dt,i)

mx=[10 10 15];
my=[10 -10 0];
mz=[10 20 -10];
Beacons=[mx' my' mz'];
mew=mew0;

for j=1:i
dt=Dt(j);
ut=U(j,:);
mew = mew + ut'*dt;
end

for j=1:3
q(j)=norm([mx(j)-mew(1) my(j)-mew(2) mz(j)-mew(3)]);
end

z=q;

end