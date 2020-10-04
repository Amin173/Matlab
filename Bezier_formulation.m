clc
clear
syms t P0x P1x P2x P3x P0y P1y P2y P3y
P0= [P0x;P0y];
P1= [P1x;P1y];
P2= [P2x;P2y];
P3= [P3x;P3y];
% B = simplify((1-t)^3 * P0 + 3*t*(1 - t^2)*P1 + 3*(1-t)*t^2*P2 + t^3*P3);
B = (1-t) * ((1-t)*P0 + t*P1) + t*((1 - t)*P1 + t*P2);
dB = simplify(diff(B, t));
f = simplify((dB(1)^2 + dB(2)^2)^0.5);
DfDP0x= simplify(diff(f,P0x));
DfDP0y= simplify(diff(f,P0y));
DfDP1x= simplify(diff(f,P1x));
DfDP1y= simplify(diff(f,P1y));
DfDP2x= simplify(diff(f,P2x));
DfDP2y= simplify(diff(f,P2y));
DfDP3x= simplify(diff(f,P3x));
DfDP3y= simplify(diff(f,P3y));

disp(DfDP0x)
disp(DfDP0y)
disp(DfDP1x)
disp(DfDP1y)
disp(DfDP2x)
disp(DfDP2y)
disp(DfDP3x)
disp(DfDP3y)