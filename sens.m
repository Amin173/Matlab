function [ y ] = sens(world,x,sensorRight,measurement)

y=x.*(world==measurement)*sensorRight+x.*(world~=measurement)*(1-sensorRight);
disp(y)
end

