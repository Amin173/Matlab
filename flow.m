clc
clear
alpha=[0.3;0.09;0.12;0.3;0.3;0.3;0.09;0.3;0.12;0.3;0.06;0.3;0.09;0.06;0.3;0.3;0.09;0.09];
p_bar=[100;50;50;100;50;50;100;50;50];

%C1 Matrix
c1=zeros(18,18);
c1(1,2)=-1;c1(1,1)=1;
c1(2,3)=-1;c1(2,2)=1;
c1(3,4)=-1;c1(3,3)=1;
c1(4,5)=-1;c1(4,3)=1;
c1(5,6)=-1;c1(5,5)=1;
c1(6,8)=-1;c1(6,7)=1;
c1(7,9)=-1;c1(7,8)=1;
c1(8,10)=-1;c1(8,9)=1;
c1(9,11)=-1;c1(9,9)=1;
c1(10,12)=-1;c1(10,11)=1;
c1(11,11)=-1;c1(11,5)=1;
c1(12,14)=-1;c1(12,13)=1;
c1(13,15)=-1;c1(13,14)=1;
c1(14,16)=-1;c1(14,14)=1;
c1(15,17)=-1;c1(15,15)=1;
c1(16,18)=-1;c1(16,16)=1;
c1(17,2)=-1;c1(17,15)=1;
c1(18,8)=-1;c1(18,16)=1;

%C2 matrix
c2=zeros(9,18);
joints=[2 1 17;
 13 15 17;
 12 13 14;
 14 18 16;
 7 18 6;
 7 9 8;
 2 3 4;
 4 5 11;
 10 11 9];
 for i=1:size(joints,1)
c2(i,[joints(i,:)])=[-1 1 1];
 end


%C3 Matrix
c3=zeros(9,18);
nodes=[1 4 6 7 10 12 13 17 18];
for i=1:length(nodes)
    c3(i,nodes(i))=1;
end

A=diag(alpha);

M=[c1 -A;
   zeros(size(c2,1),size(c1,2)) c2;
   c3 zeros(size(c3,1),size(c2,2))];

b=[zeros(size(M,1)-length(p_bar),1);p_bar];
soln=M\b;
P=soln(1:18);
Flow=soln(19:36);
Node=[1:18]';
Segment=[1:18]';
% SolutionTable=table(Node,P,Segment,Flow);
% disp(SolutionTable);

%Using reduced size system
c4=(c2/A)*c1;
M2=[c4;c3];
b2=[zeros(size(M2,1)-length(p_bar),1);p_bar];
soln_reduced=M2\b2;
P_reducedSys=soln_reduced(1:18);
Flow_reducedSys=(A\c1)*P_reducedSys;
SolutionTable2=table(Node,P,P_reducedSys,Segment,Flow,Flow_reducedSys);
disp(SolutionTable2);

%Computing the residual norm
disp('Residual vector infinity norms:')
disp('Complete solution:')
r=M*soln-b;
disp(norm(r,inf))
%Computing the residual norm of the reduced system
disp('Reduced system solution:')
r2=M2*P_reducedSys-b2;
disp(norm(r2,inf))


%% Simulatioin of Occlusion 
disp('*****Simulation of Occlusion*****');
alpha(11)=1e15;
clear r r2 soln soln_reduced M M2 b b2 c4 A
A=diag(alpha);

M=[c1 -A;
   zeros(size(c2,1),size(c1,2)) c2;
   c3 zeros(size(c3,1),size(c2,2))];

b=[zeros(size(M,1)-length(p_bar),1);p_bar];
soln=M\b;
P=soln(1:18);
Flow=soln(19:36);
Node=[1:18]';
Segment=[1:18]';
% SolutionTable=table(Node,P,Segment,Flow);
% disp(SolutionTable);


%Using reduced size system
c4=(c2/A)*c1;
M2=[c4;c3];
b2=[zeros(size(M2,1)-length(p_bar),1);p_bar];
soln_reduced=M2\b2;
P_reducedSys=soln_reduced(1:18);
Flow_reducedSys=(A\c1)*P_reducedSys;
SolutionTable2=table(Node,P,P_reducedSys,Segment,Flow,Flow_reducedSys);
disp(SolutionTable2);

%Computing the residual norm 
disp('Residual vector infinity norms:')
disp('Complete solution:')
r=M*soln-b;
disp(norm(r,inf))
%Computing the residual norm of the reduced system
disp('Reduced system solution:')
r2=M2*P_reducedSys-b2;
disp(norm(r2,inf))
%% Small resistance
disp('*****Small Resistance*****');
alpha(11)=0;
clear r r2 soln soln_reduced M M2 b b2 c4 A
A=diag(alpha);

M=[c1 -A;
   zeros(size(c2,1),size(c1,2)) c2;
   c3 zeros(size(c3,1),size(c2,2))];

b=[zeros(size(M,1)-length(p_bar),1);p_bar];
soln=M\b;
P=soln(1:18);
Flow=soln(19:36);
Node=[1:18]';
Segment=[1:18]';
% SolutionTable=table(Node,P,Segment,Flow);
% disp(SolutionTable);

%Using reduced size system
c4=(c2/A)*c1;
M2=[c4;c3];
b2=[zeros(size(M2,1)-length(p_bar),1);p_bar];
soln_reduced=M2\b2;
P_reducedSys=soln_reduced(1:18);
Flow_reducedSys=(A\c1)*P_reducedSys;
SolutionTable2=table(Node,P,P_reducedSys,Segment,Flow,Flow_reducedSys);
disp(SolutionTable2);

%Computing the residual norm
disp('Residual vector infinity norms:')
disp('Complete solution:')
r=M*soln-b;
disp(norm(r,inf));
%Computing the residual norm of the reduced system
disp('Reduced system solution:')
r2=M2*P_reducedSys-b2;
disp(norm(r2,inf));