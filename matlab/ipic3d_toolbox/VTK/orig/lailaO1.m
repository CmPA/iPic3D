addpath('/home/gianni/matlab2/matlab-parsek');
global contours

contours = 1;

[V,Bx,By,Bz,dx,dy,dz]=read_vtk('/data/gianni/paraview/lailaO1/B_xyz_cycle13000.vtk');

[nx,ny]=size(Bx);
[xx yy]=meshgrid(1:ny,1:nx);

[V,Ex,Ey,Ez,dx,dy,dz]=read_vtk('/data/gianni/paraview/tred54/E_AVG_xy_cycle13000.vtk');

% N.B. x and y are exhanged, damned matlab

lmt=round([nx/8 nx-nx/8 1 ny]);


xx=xx*dx;
yy=yy*dy;

ay=vecpot(xx,yy,Bx,By);


%disp('f1')
%h=figure(1);
%set(h,'Position' , [5 5 560 420]);
%coplot(xx,yy,EJe+dEJe+EJi+dEJi,ay,'x/d_i','y/d_i','EJTOT',lmt,[-2e-8 2e-8])
%caxis([-2 2]*1e-8)
%close all



