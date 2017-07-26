addpath '/home/gianni/Documents/matlab/matlab3'
clear all
close all


for cycle=000:10000:65000

ncycle=num2str(cycle)
ncycle1=num2str(cycle,'%06d')

% for HRmaha3D1:
time=60*(cycle/75000.0) *4 %times four to correct for change in dt between 2D and 3D
% for HRmaha3D1.v2:
% time=60*(cycle/75000.0) *2 %times two to correct for change in dt between 2D and 3D
ntime=num2str(time,'%5.2f')

dir='/data2/gianni/gda/HRmaha3D2/'


global blowup contours
blowup=0;
contours=1;

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
Lz=double(hdf5read(filename,'/collective/Lz'));
Nx=double(hdf5read(filename,'/collective/Nxc')); 
Ny=double(hdf5read(filename,'/collective/Nyc'));
Nz=double(hdf5read(filename,'/collective/Nzc'));
dx=Lx/Nx;
dy=Ly/Ny;
dz=Lz/Nz;
%xc=linspace(-45, -15, Nx);
%yc=linspace(-9, 3, Ny);
xc=linspace(0, Lx, Nx);
yc=linspace(0, Ly, Ny);
zc=linspace(0, Lz, Nz);
x=[-15 -45]; 
y=[-9 3];

for iz=1:round(Nz/7):Nz


close all
file=[dir 'B_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Bx=fread(fid,'real*8');
fclose(fid);
Bx=reshape(Bx,Nx,Ny,Nz);
Bx=squeeze(Bx(:,:,iz));

file=[dir 'B_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
By=fread(fid,'real*8');
fclose(fid);
By=reshape(By,Nx,Ny,Nz);
By=squeeze(By(:,:,iz));

file=[dir 'B_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Bz=fread(fid,'real*8');
fclose(fid);
Bz=reshape(Bz,Nx,Ny,Nz);
Bz=squeeze(Bz(:,:,iz));



titolo =[ 'time = ' ntime '   Ygsm=' num2str(gsmz2y(zc(iz)))];
nome=['Bl' ncycle1 'IZ=' num2str(iz)];
immagine_new(x,y,-Bx*code_B,nome,[-25 25],0,titolo)
nome=['Bn' ncycle1 'IZ=' num2str(iz)];
immagine_new(x,y,By*code_B,nome,[-10 10],0,titolo)
title(['time = ' ntime '   Ygsm=' num2str(gsmz2y(zc(iz)))])
end
end
