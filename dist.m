function d=dist(v)
d=(v(1,:).^2+v(2,:).^2).^0.5;
d=[d;d];
end