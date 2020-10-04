clc
clear
mu=[0 0 0];
SIGMA=[0.01 0 0;0 0.01 0;0 0 1e4];

th=0:pi/4:2*pi-pi/4;
x1=-5:5;
x2=-5:5;
n=length(x1);
X1=x1'*ones(1,n);
X2=meshgrid(x2);

X1=reshape(X1,[n^2,1]);
X2=reshape(X2,[n^2,1]);


%generate cell centers
% xc1=0.5:1:4.5;
% xc2=0.5:1:4.5;
% m=length(xm1);
% Xc1=xc1'*ones(1,m);
% Xc2=meshgrid(xc2);
% 
% Xc1=reshape(Xc1,[m^2,1]);
% Xc2=reshape(Xc2,[m^2,1]);

for i=1:length(th)
    th_k=th(i)*ones(n^2,1);
X=[X1 X2 th_k];

p{i}= mvncdf(X,mu,SIGMA);
end

for j=1:length(th)
p_j=reshape(p{j},[n,n]);
pd(:,1,j)=p_j(:,1);
for i2=2:n
pd(:,i2,j)=p_j(:,i2)-p_j(:,i2-1);
end
pd2(1,:,j)=pd(1,:,j);
for i1=2:n
    pd2(i1,:,j)=pd(i1,:,j)-pd(i1-1,:,j);
end
end
pd2=pd2/sum(sum(sum(pd2)));  %initial probability mass distribution

p_bar=pd2;
pd=zeros(n,n,length(th));
pd_ij=zeros(n,n);
for k=1:length(th)
    for i=1:n
        for j=1:n
            x=x1(i)+cos(th(k))
            y=x2(j)+sin(th(k))
            dist=sqrt((X1-x).^2+(X2-y).^2);
            
            i1=(x1==X1(dist==min(dist)))
            i2=(x2==X2(dist==min(dist)))
            pd_ij([i1],[i2])=p_bar(i,j,k)/(sum(i1)*sum(i2));
            pd(:,:,k)=pd(:,:,k)+pd_ij;
        end
    end
end
         pd=pd/sum(sum(sum(pd)));