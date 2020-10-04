function [x_extreme] = find_extreme2(dy)
j=zeros(length(dy),1);
for i=2:length(x)-1
    j(i)=dy(i)*dy(i+1)<0;
end
j(1)=dy(1)==0;
j(end)=dy(end)==0;
x_extreme=j==1;
x_extreme=reshape(x_extreme,[length(x_extreme),1]);

end