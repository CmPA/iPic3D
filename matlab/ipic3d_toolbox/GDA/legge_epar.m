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
Lx=hdf5read(filename,'/collective/Lx'); 
Ly=hdf5read(filename,'/collective/Ly');
Nx=hdf5read(filename,'/collective/Nxc'); 
Ny=hdf5read(filename,'/collective/Nyc');
xc=linspace(0, Lx, Nx);
yc=linspace(0, Ly, Ny);
x=[-15 -45]; 
y=[-9 3];

close all
file=[dir 'E_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Ex=fread(fid,'real*8');
fclose(fid);
Ex=reshape(Ex,Nx,Ny);
file=[dir 'E_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Ey=fread(fid,'real*8');
fclose(fid);
Ey=reshape(Ey,Nx,Ny);
file=[dir 'E_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Ez=fread(fid,'real*8');
fclose(fid);
Ez=reshape(Ez,Nx,Ny);

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
Epar=(Ex.*Bx+Ey.*By+Ez.*Bz)./sqrt(B2);
Epar=Epar(10:end-10,10:end-10);
Ay=vecpot_uniform(xc,yc,Bx,By);
immagine(x,y,Epar*code_E,['Epar' ncycle],[min(Epar(:)) max(Epar(:))]*code_E,0)
print('-depsc','-r300','epar.eps')
hold on
contour(-15-xc/Lx*30,i-9+yc/Ly*12,Ay',50,'w')
set(gcf, 'Renderer', 'zbuffer');
%ylim([-1 0.7])
%xlim([-21.5 -17.5])
ylim([1 2])
xlim([-15.5 -15])
print('-depsc','-r300','Eparcombo.eps')
print('-dpng','-r300','Eparcombo.png')



