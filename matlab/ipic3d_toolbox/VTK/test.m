addpath('/home/gianni/matlab2/matlab-parsek');
global contours

contours = 1;

for it=19000:1000:19000

[V,Jx,Jy,Jz,dx,dy,dz]=read_vtk(['/shared02/gianni/maha2/data1/Je_xyz_cycle' num2str(it) '.vtk']); 
[nx,ny]=size(Jx);
[xx yy]=meshgrid(1:ny,1:nx);

% N.B. x and y are exhanged, damned matlab

lmt=round([nx/8 nx-nx/8 1 ny]);


%
% Read electron pressure tensor
%

xx=xx*dx;
yy=yy*dy;

%
% Pressure tensor
%

ay=Jex*0;
disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(it,xx,yy,Jex',ay,'x/d_i','y/d_i','Jex',lmt)
close all


end




