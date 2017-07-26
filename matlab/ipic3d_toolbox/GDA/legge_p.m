clear all
close all

dir='/shared02/gianni/maha2/data4/'
cycle=75000
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

legge_P=1


if(legge_P)
close all
file=[dir 'Pe_per1_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Peper1=fread(fid,'real*8');
fclose(fid);
Peper1=reshape(Peper1,Nx,Ny);
file=[dir 'Pe_per2_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Peper2=fread(fid,'real*8');
fclose(fid);
Peper2=reshape(Peper2,Nx,Ny);

file=[dir 'rho_0_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Ne=fread(fid,'real*8');
fclose(fid);
Ne=reshape(Ne,Nx,Ny);

Teper=(Peper1+Peper2)./(-Ne)/2;


%h = fspecial('gaussian',[5 5], .5);
%J = filter2(h, Je+Ji);

TeoTi=1/5;
immagine(x,y,Teper*code_T/e/1e3/TeoTi,['Teper' ncycle],[1 8]/TeoTi/2*1.5,3)
%immagine(x,y,Teper*9/1.977e-6,'Teper',[1 30],3)



end
