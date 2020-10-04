clc
clear
rp = 0.005;
rt = 0.04;
L = 2* 6* 0.0254;
Vt= pi * rt^2 * L;
Vp= 4/3 * pi * rp^3;
m = Vp/Vt;
%%
figure(1)
hold on
for n=1:500
    plot(n, 100*(n* m ) , '.')
end
xlabel('Number of particles')
ylabel('Particle volume fraction %')
grid on
grid minor 
%%
figure(2)
hold on
for i = 20:10:40
    plot(i, round(i/(100*m)),'s', 'LineWidth', 3)
end
xlabel('Particle volume fraction %')
ylabel('Nummber of particles')
grid on
grid minor
