clear all
close all
global Ex Ey Ez Bx By Bz xc yc qom


cycle=37000
ncycle=num2str(cycle)
ncycle1=num2str(cycle,'%06d')

dir='/shared02/gianni/maha2/data1/'
if(cycle>19000)
dir='/shared02/gianni/maha2/data2/'
end
if (cycle>37000)
dir='/shared02/gianni/maha2/data3/'
end
if (cycle>55400)
dir='/shared02/gianni/maha2/data4/'
end

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


qom = -256;


addpath '/home/gianni/matlab2/matlab-parsek'
filename=[dir 'settings.hdf'];
Lx=double(hdf5read(filename,'/collective/Lx')); 
Ly=double(hdf5read(filename,'/collective/Ly'));
Nx=double(hdf5read(filename,'/collective/Nxc')); 
Ny=double(hdf5read(filename,'/collective/Nyc'));
dx=Lx/Nx;
dy=Ly/Ny;
%xc=linspace(-45, -15, Nx);
%yc=linspace(-9, 3, Ny);
xc=linspace(0, Lx, Nx);
yc=linspace(0, Ly, Ny);
x=[-15 -45]; 
y=[-9 3];


close all
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

Ay=zeros(size(Bx));
Ay=vecpot_uniform(xc,yc,Bx*dy/dx,By);

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

Nsm=0;
Ex=smooth(Ex,Nsm);
Ey=smooth(Ey,Nsm);
Ez=smooth(Ez,Nsm);
Bx=smooth(Bx,Nsm);
By=smooth(By,Nsm);
Bz=smooth(Bz,Nsm);

%options=odeset('AbsTol',10e-1,'RelTol',10e-1,'Stats','off'); 
counter=0
for vxp=-.1:.01:.1
[t,y]=ode45(@newton,0:.1:2000 ,[24 15 0 vxp .01 .01])
%[t,y]=ode45(@newton,0:.1:2000 ,[35 15 0 vxp .01 .01])
counter=counter+1
save(['autotraject' num2str(counter) '.mat'],'t','y')
%plot(y(:,1),y(:,2))
plot_traj
end
