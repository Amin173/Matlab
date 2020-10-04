clc
clear
n=10;
x=1:n;

a=ones(n,1);
for j=1:n
    for k=1:n
        if k~=j
        a(j)=a(j)*(x(j)-x(k));
        end
    end
end

D=zeros(n,n);
for i=1:9
    for j=1:9
        if j~=i
        D(i,j)=a(i)/(a(j)*(x(i)-x(j)));
        elseif j==i
            for k=1:n
                if k~=j
            D(j,j)=D(j,j)+(1/(x(j)-x(k)));
                end
            end
        end
    end
end
%%
xx=3;
px=ones(n,1);
for j=1:n
        for k=1:n
        if k~=j
    px(j)=(1/a(j))*(xx-x(k));
        end
        end
end

dpx=zeros(n,1);
for j=1:n
        for k=1:n
        if k~=j
    dpx(j)=px(j)*(xx-x(k))+dpx(j);
        end
        end
end