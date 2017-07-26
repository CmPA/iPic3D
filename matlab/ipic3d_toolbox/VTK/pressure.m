addpath('/home/gianni/matlab2/matlab-parsek');
global contours

contours = 1;

for it=0:1000:20000

[V,Bx,By,Bz,dx,dy,dz]=read_vtk(['/shared02/gianni/tred60/data/B_AVG_xy_cycle' num2str(it) '.vtk']);

[nx,ny]=size(Bx);
[xx yy]=meshgrid(1:ny,1:nx);

% N.B. x and y are exhanged, damned matlab

lmt=round([nx/8 nx-nx/8 1 ny]);


%
% Read electron pressure tensor
%

xx=xx*dx;
yy=yy*dy;

ay=vecpot(xx,yy,Bx,By);

[Pexx,Pexy,Pexz,Peyy,Peyz,Pezz,TPAR, TPER1, TPER2, EPS]=read_vtk_multiscalar(['/data/gianni/paraview/tred60/AVGPe_AVG_xy_cycle' num2str(it) '.vtk'],10);
[Rhoe,Rhoi]=read_vtk_multiscalar(['/shared02/gianni/tred60/data/AVGrho_AVG_xy_cycle' num2str(it) '.vtk'],2);
Rhoe=-Rhoe;


load('cm_new')

%
% Pressure tensor
%

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(it,xx,yy,Pexx',ay,'x/d_i','y/d_i','Pexx',lmt)
close all

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(it,xx,yy,Peyy',ay,'x/d_i','y/d_i','Peyy',lmt)
close all

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(it,xx,yy,Pezz',ay,'x/d_i','y/d_i','Pezz',lmt)
close all


disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(it,xx,yy,Pexy',ay,'x/d_i','y/d_i','Pexy',lmt)
close all

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(it,xx,yy,Pexz',ay,'x/d_i','y/d_i','Pexz',lmt)
close all

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(it,xx,yy,Peyz',ay,'x/d_i','y/d_i','Peyz',lmt)
close all




disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(it,xx,yy,TPAR'./Rhoe',ay,'x/d_i','y/d_i','Tpar',lmt)
close all

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(it,xx,yy,TPER1'./Rhoe',ay,'x/d_i','y/d_i','Tperp1',lmt)
close all

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(it,xx,yy,TPER2'./Rhoe',ay,'x/d_i','y/d_i','Tperp2',lmt)
close all

disp('f1')
h=figure(1);
set(h,'Position' , [5 5 560 420]);
coplot(it,xx,yy,EPS',ay,'x/d_i','y/d_i','EPS',lmt)
close all


end




