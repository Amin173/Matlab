clc; clear;
t = zeros(1, 12);
t(1:6) = pi;

A = zeros(2);
B = [cos(t) ; sin(t)];
Q = 0.0001* eye(2);

R = 2 * eye(12);
R(1, [3, end-1]) = -1;
R(end, [end-2, 2]) = -1;
R(2, [4, end]) = -1;
R(end-1, [end-3, 1]) = -1;
for i= 3:length(t)-2
    R(i, i-2) = -1;
    R(i, i+2) = -1;
end

[K,S,e] = lqrd(A,B,Q,R,1)