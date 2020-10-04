l = 102*10^(-3); T = 100; D = 10*10^(-3); r= D/2; row = 1000; lan = 2*10^6;miu = 5000; f = 1000; 
E = 15000; nu = 0.5; G = 5000; 
omega = 2*pi*f;
Area = pi*r^2; In = pi/4*r^4;
dl = T*l/(Area*E);
c = l+dl;
dD = -nu*D*dl/l;
d = D-dD;
longitudinal = sqrt(E/row);
longitudinal_damping = sqrt(15+3j);
torsional = sqrt(G/row);
string = sqrt(T/(Area*row));
thin_beam = sqrt(omega)*(E*In/(row*Area))^(1/4);
% Thick beam
kapa = 10/9;
a = (row*omega^2*(G+(kapa*E)))/(G*E);
b = (4*(row^2)*(omega^2)*kapa*((omega^2)-(G*Area/(row*In*kapa))))/(E*G);
k = -1*sqrt((-a-sqrt((a^2)-b))/2);
thick_beam = omega/(k);
%pre-stress 
alfa = sqrt(((-2*T)/(E*pi*r^4))+(sqrt(((2*T/(E*pi*r^4))^2)+((4*row*(omega)^2)/(E*r^2)))));
pre_stress = omega/alfa;
