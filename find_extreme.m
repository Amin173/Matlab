function [x_extreme,y_extreme] = find_extreme(x,y)
j=zeros(length(x),1);
for i=2:length(x)-1
    j(i)=((y(i)>y(i-1))&&(y(i+1)<y(i)))||((y(i)<y(i-1))&&(y(i+1)>y(i)));
end
x_extreme=x(j==1);
y_extreme=y(j==1);
x_extreme=reshape(x_extreme,[length(x_extreme),1]);
y_extreme=reshape(y_extreme,[length(x_extreme),1]);

end

