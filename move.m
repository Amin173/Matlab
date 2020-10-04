function [ y ] = move(n,x,m,pMove)

for i=1:n
    for j=1:n
        i1=mod(i+m(1),n);
        i1=i1*(i1>0)+(i1+n)*(i1<=0);
        i2=mod(j-m(2),n);
        i2=i2*(i2>0)+(i2+n)*(i2<=0); 
        y(i,j)=pMove*x(i1,i2)+(1-pMove)*x(i,j);
    end
end
disp(y)
end

