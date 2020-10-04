
clc; clear; close all;
for n=1:5
N=10; %number of nodes (robots)
tpn=n*ones(N,1); %tendons per node
c=zeros(N); %connecivity matrix 
% If c(i,j)=1 node i is connected to node j 

for i=1:N
    for j=1:tpn(i)
        if (i+j)<=N
            c(i+j,i)=1;
        else
            c(j+i-N,i)=1;
        end        
    end
end

f=zeros(N,N); %force matrix
% Columns represent the servo moounted robot and rows reprent the robots that the servo is cabled to

for i=1:N
    
    fc=1/tpn(i); %foce per cable
    NN=[1:N 1:N 1];
        for  j=i+tpn(i):-1:i+1
         f(NN(j),NN(i))=fc*(c(NN(j),NN(i))==1)+f(NN(j+1),NN(i));
        end
        
    end
    total_nodal_force=sum(f);  % The force at each node in [unit servo force]
    normalized_nodal_force=sum(f)/N;
    
%     disp('The force at each node in [unit servo force]:')
%     disp(total_nodal_force)

tnf(n)=total_nodal_force(1);
end

plot(1:5,tnf,'+:','LineWidth',2)
xlabel('Number of cables per servo')
ylabel('Totoal force in each segment [unit servo force]')
set(gca,'xtick',1:5)