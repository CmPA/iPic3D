clear all
close all

dir='/shared02/gianni/maha2/data2/'
cycle=37000
ncycle=num2str(cycle);

global blowup contours
blowup=0;
contours=0;

code_E = 2060.21;
code_B = 6.87213e-06;
code_B=code_B *1e9; % to convert from Tesla to nT
code_J = 1.20082e-05;
code_J = code_J*1e9; % to convert to nA/m^2
code_V = 2.99792e+08;
code_V=code_V/1e3; %to convert to Km/s
code_T = 1.50326e-10;
code_n = 0.25;
e=1.6e-19;


addpath '/home/gianni/matlab2/matlab-parsek'
filename=[dir 'settings.hdf'];
Lx=double(hdf5read(filename,'/collective/Lx'));
Ly=double(hdf5read(filename,'/collective/Ly'));
Nx=double(hdf5read(filename,'/collective/Nxc'));
Ny=double(hdf5read(filename,'/collective/Nyc'));
dx=Lx/Nx;
dy=Ly/Ny;
xc=linspace(0, Lx, Nx);
yc=linspace(0, Ly, Ny);
x=[-15 -45]; 
y=[-9 3];

close all
file=[dir 'Je_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Ex=fread(fid,'real*8');
fclose(fid);
Ex=reshape(Ex,Nx,Ny);
file=[dir 'Je_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Ey=fread(fid,'real*8');
fclose(fid);
Ey=reshape(Ey,Nx,Ny);
file=[dir 'Je_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Ez=fread(fid,'real*8');
fclose(fid);
Ez=reshape(Ez,Nx,Ny);
file=[dir 'rho_0_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
rhoe=fread(fid,'real*8');
fclose(fid);
rhoe=reshape(rhoe,Nx,Ny);
file=[dir 'rho_1_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
rhoi=fread(fid,'real*8');
fclose(fid);
rhoi=reshape(rhoi,Nx,Ny);

file=[dir 'B_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Bx=fread(fid,'real*8');
fclose(fid);
Bx=reshape(Bx,Nx,Ny);
file=[dir 'B_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
By=fread(fid,'real*8');
fclose(fid);
By=reshape(By,Nx,Ny);
file=[dir 'B_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Bz=fread(fid,'real*8');
fclose(fid);
Bz=reshape(Bz,Nx,Ny);

B2=(Bx.*Bx+By.*By+Bz.*Bz);
va=sqrt(B2./4./pi./rhoi);
Vex=Ex./rhoe./va;
Vey=Ey./rhoe./va;
Vez=Ez./rhoe./va;

Ay=vecpot_uniform(xc,yc,Bx*dy/dx,By);

immagine(x,y,Vex,['VexOVa' ncycle],[-5 5],5)
%print('-depsc','-r300',['vexova' ncycle '.eps'])
hold on
contour(-15-xc/Lx*30,i-9+yc/Ly*12,Ay',50,'w')
set(gcf, 'Renderer', 'zbuffer');
%ylim([-5 -1])
%xlim([-35 -20])
print('-depsc','-r300',['vexovacombo' ncycle '.eps'])
print('-dpng','-r300',['vexovacombo' ncycle '.png'])

immagine(x,y,Vey,['VeyOVa' ncycle],[-5 5],5)
%print('-depsc','-r300',['veyova' ncycle '.eps'])
hold on
contour(-15-xc/Lx*30,i-9+yc/Ly*12,Ay',50,'w')
set(gcf, 'Renderer', 'zbuffer');
%ylim([-5 -1])
%xlim([-35 -20])
print('-depsc','-r300',['veyovacombo' ncycle '.eps'])
print('-dpng','-r300',['veyovacombo' ncycle '.png'])


immagine(x,y,Vez,['VezOVa' ncycle],[-5 5],5)
%print('-depsc','-r300',['vezova' ncycle '.eps'])
hold on
contour(-15-xc/Lx*30,i-9+yc/Ly*12,Ay',50,'w')
set(gcf, 'Renderer', 'zbuffer');
%ylim([-5 -1])
%xlim([-35 -20])
print('-depsc','-r300',['vezovacombo' ncycle '.eps'])
print('-dpng','-r300',['vezovacombo' ncycle '.png'])

return
% below a failed attempt to plot the velocity field. Too noisy and too structured.
close all
Nsm=20
Vex=smooth(Ex,Nsm);
Vey=smooth(Ey,Nsm);
hlines=streamslice(xc,yc,-fliplr(Vex')*dy/dx,fliplr(Vey'));
set(hlines,'Color','r')
set(gcf, 'Renderer', 'zbuffer');

print('-dpng','-r300',['vestream' ncycle '.png'])


