clc
clear
syms tfinal Vmax V0 tmax d K
Vfinal = V0 + K*tfinal;
dt1= tfinal;
dt2 = d/Vfinal;
dt(tfinal) = dt1 + dt2;
dt_opt = simplify(solve(diff(dt, tfinal)==0, tfinal));
simplify(dt_opt(2))
%%
clear
clc
syms dt m n dd V0 K V
T(m) = m*dt + (n-m)*dd/(K*m*dt+V0);
mopt(dt, dd)= solve(diff(T)==0,m);
simplify(mopt(1, V));