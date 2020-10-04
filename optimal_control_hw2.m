clc
clear
%% Problem 2
E = [1 0 1 0 0; 0 0 1 0 -.0423;0 0 -.106 0 1; 
    0 0 0 1 0; 0 1 0 0 0];
A = [-.0297 0 0 .0438 0;
    .379 0 -.0096 0 -.0125;
    -1.17 0 .129 0 -.79;
    0 0 0 0 1;
    0 0 1 0 0];
dr = 1;
da = 1;
phi0 = 1;
e0 =1;
u = [dr;da];
B = [0 0;
    1 0;
    0 1;
    0 0;
    0 0];
R = [1/dr^2 0;
    0 1/da^2];
Q = [1/e0^2 1 0 0 0;1 1/e0^2 0 0 0; 0 0 0 0 0; 0 0 0 1/phi0^2 0; 0 0 0 0 0];

[K,S,e] = lqr(inv(E)*A, inv(E)*B, Q, R, 0);

x0 = [0;0;0.1;0;0.2]; % initial state
dt = 0.1;
n = 100;
% for i = 1:n
% u = -K * x0;
% xdot = inv(E) * A * x0 + inv(E) * B * u;
% x(:,i) = x0 + xdot*dt;
% x0 = x(:,i);
% end

G = ss(inv(E) * A,inv(E) * B,eye(length(x0)), 0*B);
sys = feedback(G,K);
initial(sys,x0)
% plot((1:n)/dt, x(:,:)')
% xlabel('Time [s]')
% ylabel('System response')
% legend('x1:beta', 'x2: sai', 'x3: r', 'x4: phi', 'x5: p')
%% Problem 3 (a)
clc; clear;
R = 1;
q11 = 1;
Q = [q11 0 0; 0 0 0; 0 0 0];
A = [ 0 1 0; 8.75e4 0 7e4; 0 -0.5/0.4 -2.2/0.00368];
B = [0;0;1];
[K,S,e] = lqr(A, B, Q, R, 0);
G = ss(A, B, eye(3),zeros(3,1));
sys = feedback(G,K);
pole(sys)
%% (b) 
q11 = [0.001, .01, .1, 1, 10, 100, 1000];
for i=1:length(q11)
    x0 = [0.1;0;0];
    R = 1;
    Q = [q11(i) 0 0; 0 0 0; 0 0 0];
    A = [ 0 1 0; 8.75e4 0 7e4; 0 -0.5/0.4 -2.2/0.00368];
    B = [0;0;1];
    [K,S,e] = lqr(A, B, Q, R, 0);
    G = ss(A, B, eye(3),zeros(3,1));
    sys = feedback(G,K);
    res = initial(sys,x0);
    hold on
    for j= 1:3
        subplot(4,1,j)
        plot(res(:,j))
    end
    subplot(4,1,4)
    plot(-0.00368*K*res')
    rms_vc(i) = rms(-K*res');
    rms_x(i) = rms(res(:,1));
end
xlabel('Time [s]')
subplot(4,1,1)
ylabel('x[mm]')
subplot(4,1,2)
ylabel('x_{dot}[mm/s]')
subplot(4,1,3)
ylabel('Ic[A]')
subplot(4,1,4)
ylabel('Vc[V]')
figure()
loglog(q11,rms_vc, q11, rms_x)
xlabel('q11')
legend('RMS(u)','RMS(x)')
%% (c) 
q11 = 1;
q22 = [0.001, .01, .1, 1, 10, 100, 1000];
for i=1:length(q11)
    x0 = [0.1;0;0];
    R = 1;
    Q = [q11 0 0; 0 q22(i) 0; 0 0 0];
    A = [ 0 1 0; 8.75e4 0 7e4; 0 -0.5/0.4 -2.2/0.00368];
    B = [0;0;1];
    [K,S,e] = lqr(A, B, Q, R, 0);
    G = ss(A, B, eye(3),zeros(3,1));
    sys = feedback(G,K);
    res = initial(sys,x0);
    hold on
    for j= 1:3
        subplot(4,1,j)
        plot(res(:,j))
    end
    subplot(4,1,4)
    plot(-0.00368*K*res')
    rms_vc(i) = rms(-K*res');
    rms_x(i) = rms(res(:,1));
end
xlabel('Time [s]')
subplot(4,1,1)
ylabel('x[mm]')
subplot(4,1,2)
ylabel('x_{dot}[mm/s]')
subplot(4,1,3)
ylabel('Ic[A]')
subplot(4,1,4)
ylabel('Vc[V]')
figure()
loglog(q22,rms_vc, q22, rms_x)
xlabel('q22')
legend('RMS(u)','RMS(x)')