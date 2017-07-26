addpath('/home/gianni/matlab2/matlab-parsek');
global contours

contours = 1;

[V,Bx,By,Bz,dx,dy,dz]=read_vtk('/data/gianni/paraview/tred60/B_AVG_xy_cycle16000.vtk');

[nx,ny]=size(Bx);
[xx yy]=meshgrid(1:ny,1:nx);

[V,Ex,Ey,Ez,dx,dy,dz]=read_vtk('/data/gianni/paraview/tred60/E_AVG_xy_cycle16000.vtk');

% N.B. x and y are exhanged, damned matlab

lmt=round([nx/8 nx-nx/8 1 ny]);


xx=xx*dx;
yy=yy*dy;

ay=vecpot(xx,yy,Bx,By);

[V,dPx,dPy,dPz,dx,dy,dz]=read_vtk('/data/gianni/paraview/tred60/AVGdelPoy_AVG_xy_cycle16000.vtk');

[V,Px,Py,Pz,dx,dy,dz]=read_vtk('/data/gianni/paraview/tred60/AVGJPoy_AVG_xy_cycle16000.vtk');

load('cm_new')

%
% Electric field
%
disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(xx,yy,Ex,ay,'x/d_i','y/d_i','Ex',lmt)
close all

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(xx,yy,Ey,ay,'x/d_i','y/d_i','Ey',lmt)
close all

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(xx,yy,Ez,ay,'x/d_i','y/d_i','Ez',lmt)
close all
%
% Poynting flux
%

%  fluctuating

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(xx,yy,dPx,ay,'x/d_i','y/d_i','dPoy_x',lmt)
close all

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(xx,yy,dPy,ay,'x/d_i','y/d_i','dPoy_y',lmt)
close all

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(xx,yy,dPz,ay,'x/d_i','y/d_i','dPoy_z',lmt)
close all

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(xx,yy,sqrt(dPx.^2+dPy.^2+dPz.^2),ay,'x/d_i','y/d_i','dPoy_m',lmt)
close all

%  AVERAGE

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(xx,yy,Px,ay,'x/d_i','y/d_i','Poy_x',lmt)
close all

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(xx,yy,Py,ay,'x/d_i','y/d_i','Poy_y',lmt)
close all

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(xx,yy,Pz,ay,'x/d_i','y/d_i','Poy_z',lmt)
close all

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(xx,yy,sqrt(Px.^2+Py.^2+Pz.^2),ay,'x/d_i','y/d_i','Poy_m',lmt)
close all

%
% E.J Electrons
%

%  fluctuating

[V,dPx,dPy,dPz,dx,dy,dz]=read_vtk('/data/gianni/paraview/tred60/AVGdelJdotEe_AVG_xy_cycle16000.vtk');

[V,Px,Py,Pz,dx,dy,dz]=read_vtk('/data/gianni/paraview/tred60/AVGJdotEe_AVG_xy_cycle16000.vtk');

dEJe= dPx+dPy+dPz;
EJe= Px+Py+Pz;

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(xx,yy,dEJe,ay,'x/d_i','y/d_i','dEJe',lmt,[-2e-8 2e-8])
close all

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(xx,yy,EJe,ay,'x/d_i','y/d_i','EJe',lmt,[-2e-8 2e-8])
close all

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(xx,yy,EJe+dEJe,ay,'x/d_i','y/d_i','EJeTOT',lmt,[-2e-8 2e-8])
close all


%
% E.J Ions
%

%  fluctuating

[V,dPx,dPy,dPz,dx,dy,dz]=read_vtk('/data/gianni/paraview/tred60/AVGdelJdotEi_AVG_xy_cycle16000.vtk');

[V,Px,Py,Pz,dx,dy,dz]=read_vtk('/data/gianni/paraview/tred60/AVGJdotEi_AVG_xy_cycle16000.vtk');

dEJi= dPx+dPy+dPz;
EJi= Px+Py+Pz;

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(xx,yy,dEJi,ay,'x/d_i','y/d_i','dEJi',lmt,[-2e-8 2e-8])
caxis([-2 2]*1e-8)
close all

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(xx,yy,EJi,ay,'x/d_i','y/d_i','EJi',lmt,[-2e-8 2e-8])
caxis([-2 2]*1e-8)
close all

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(xx,yy,EJi+dEJi,ay,'x/d_i','y/d_i','EJiTOT',lmt,[-2e-8 2e-8])
caxis([-2 2]*1e-8)
close all



disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(xx,yy,dEJe+dEJi,ay,'x/d_i','y/d_i','dEJ',lmt,[-2e-8 2e-8])
caxis([-2 2]*1e-8)
close all

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(xx,yy,EJe+EJi,ay,'x/d_i','y/d_i','EJ',lmt,[-2e-8 2e-8])
caxis([-2 2]*1e-8)
close all

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(xx,yy,EJe+dEJe+EJi+dEJi,ay,'x/d_i','y/d_i','EJTOT',lmt,[-2e-8 2e-8])
caxis([-2 2]*1e-8)
close all



