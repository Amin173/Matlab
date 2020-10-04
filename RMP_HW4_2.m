clc
clear
mu=[0 0 0];
SIGMA=[0.01 0 0;0 0.01 0;0 0 1e4];

th=-pi:pi/4:pi;
d=0.05;
x1=-2:d:2;
x2=-2:d:2;
n=length(x1);
% Obtaining initial probability distribution
for k=1:length(th)
[X1,X2] = meshgrid(x1,x2);
TH=th(k)*ones(length(X1(:)),1);
X = [X1(:) X2(:) TH];
y{k} = mvnpdf(X,mu,SIGMA);
y{k} = reshape(y{k},length(x2),length(x1));

end
contour3(x1,x2,y{th==0},[0.0001 0.001 0.01 0.05 0.15 0.25 0.35])
xlabel('x')
ylabel('y')
title('Initial probability distribution (theta=0)');

% Obtaining initial cumulative probability distribution over grid cells
for j=1:length(th)
for i1=1:n
for i2=1:n
    xl=[x1(i1)-d/2,x2(i2)-d/2,th(j)-1e-3];
    xu=[x1(i1)+d/2,x2(i2)+d/2,th(j)+1e-3];
    [p,err]=mvncdf([xl],[xu],mu,SIGMA);
    Q(i1,i2,j)=p;
end
end
end
Q=Q/sum(sum(sum(Q)));
figure();
contour3(x1,x2,Q(:,:,th==0));
xlabel('x')
ylabel('y')
title('Initial cumulative probability distribution (theta=0)');

Q_bar=Q;
Q_hat=zeros(n,n,length(th)); 
q_ij=zeros(n,n);  %dummy variable
for k=1:length(th)
    for i=1:n
        for j=1:n
            x=x1(i)+cos(th(k));
            y=x2(j)+sin(th(k));
            dist=sqrt((X1(:)-x).^2+(X2(:)-y).^2);
            
            i1=(x1==X1(dist==min(dist)));
            i2=(x2==X2(dist==min(dist)));
            q_ij([i1],[i2])=Q_bar(i,j,k)/(sum(i1)*sum(i2));
            Q_hat(:,:,k)=Q_hat(:,:,k)+q_ij;
        end
    end
end
         Q_hat=Q_hat/sum(sum(sum(Q_hat)));
         
         figure();
         contour3(x1,x2,Q_hat(:,:,th==pi/2));
         xlabel('x');
         ylabel('y');
         title('Cumulative probability distribution after motion (prediction step)-Theta=pi/2');

         %Measrement update
  sigma_z=0.01;
  z_t=1;
  for j=1:length(th)
  for i1=1:n
    for i2=1:n
    xl=[x1(i1)-d/2];
    xu=[x1(i1)+d/2];
    [p,err]=mvncdf([xl],[xu],z_t,sigma_z);
    Q_m(i1,i2,j)=p;
    end
  end
  end
  Q_m=Q_m/sum(sum(sum(Q_m)));
  
  %new estimate
  Q_bar=Q.*Q_m;
  Q_bar=Q_bar/sum(sum(sum(Q_bar)));
  
         figure();
         contour3(x1,x2,Q_bar(:,:,th==-pi));
         xlabel('x');
         ylabel('y');
         title('Cumulative probability distribution (correction step)-Theta=-pi');
         

%% Problem 2
clc
clear
mu=[0 0 0];
SIGMA=[0.01 0 0;0 0.01 0;0 0 1e4];

th=-pi:pi/4:pi;
d=0.05;
x1=-4:d:4;
x2=-4:d:4;
n=length(x1);

measurement=[1];
Q=0.01;

for k=1:length(th)
[X1,X2] = meshgrid(x1,x2);
TH=th(k)*ones(length(X1(:)),1);
X = [X1(:) X2(:) TH];

p= mvnpdf(X,mu,SIGMA);
x_m(:,:,k)=reshape(p,[n,n]);
end
x_m=x_m/sum(sum(sum(x_m)));
         figure();
         contour3(x1,x2,x_m(:,:,th==0));
         xlabel('x');
         ylabel('y');
         title('Initial probability distribution (theta=0)');

% prediction step
for k=1:length(th)
[X1,X2] = meshgrid(x1,x2);
TH=th(k)*ones(length(X1(:)),1);
X = [X1(:) X2(:) TH];

mu_bar=mu+[cos(th(k)) sin(th(k)) 0];
SIGMA_bar=eye(3)*SIGMA*eye(3); %nois-less motion
p= mvnpdf(X,mu_bar,SIGMA_bar); %obtaining the probability distribution over sample points
x_m(:,:,k)=reshape(p,[n,n]);
end
x_m=x_m/sum(sum(sum(x_m)));

         figure();
         contour3(x1,x2,x_m(:,:,th==pi/2));
         xlabel('y');
         ylabel('x');
         title('Probability distribution after motion (prediction)-Theta=pi/2');

X_bar=zeros(n,n,length(th));

% correction step
for i=1:length(measurement)
z=measurement(i);
w=mvnpdf(X1(:),z,Q);
w=reshape(w,n,n);
for k=1:length(th)
X_bar(:,:,k)=x_m(:,:,k).*w+X_bar(:,:,k);
end
% resampling
M=n*n*length(th);
r=rand(1)/M;
c=w(1);
i=1;
X_bar=reshape(X_bar,[M,1]);
X_t=zeros(M,1);
for m=1:M
    U=r+(m-1)/M;
    while U>c
        i=i+1;
        c=c+w(i);
    end
    X_t(m)=X_bar(i);
end
X_t=X_t/sum(sum(sum(X_t)));
X_bar=reshape(X_t,[n,n,length(th)]);
end
         figure();
         contour3(x1,x2,X_bar(:,:,th==-pi));
         colorbar;
         xlabel('x');
         ylabel('y');
         title('Highest probability after resampling- Theta=-pi')
         xlim([-4 -2])
         ylim([-4 -2])
