clc;clear;
std=20e-6*rand;
[z , PixelWidth, PSD] = artificial_surf(std, 0.1, 0.1, 512 , 512);
[n,m] = size(z);
x = linspace(0,(m-1) * PixelWidth , m);
y = linspace(0,(n-1) * PixelWidth , n);
[X,Y] = meshgrid(x,y);
mesh(X,Y,z)
colormap jet
c = colorbar;
c.Label.String = 'Height [m]';
axis equal
M=[X(:) Y(:) z(:)];
% csvwrite('rough_surface_200um_10cm_3.csv',M);
csvwrite('C:\Users\amink\Documents\MATLAB\rough_surface_200um_10cm_30.csv',M)
disp("Max height of the surface is:");
disp(max(z(:)));
%%
clc
clear

files=dir('C:\Users\amink\Documents\MATLAB\rough_surface*.csv');

for i=1:length(files)
%     eval(['load ' files(i).name  ' -ascii'])
    read=csvread(files(i).name);



% read=csvread('C:\Users\amink\Documents\MATLAB\soc-ocv\rough_surface_200um_10cm.csv');

z=read(:,3);
x=read(:,1);
y=read(:,2);

Rq=rms(z);
Ra=mean(abs(z));
Rp=max(z);
Rv=min(z);
Ry=Rp-Rv;
disp([files(i).name '    Rp: ' num2str(Rp) '    Rq: ' num2str(Rq) '   Ra: ' num2str(Ra) '   Ry: ' num2str(Ry)]); 
end
