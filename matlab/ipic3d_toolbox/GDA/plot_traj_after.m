clear all
close all
global Ex Ey Ez Bx By Bz xc yc qom Lx Ly


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
xRE=[-15 -45]; 
yRE=[-9 3];


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

Nsm=10;
Ex=smooth(Ex,Nsm);
Ey=smooth(Ey,Nsm);
Ez=smooth(Ez,Nsm);
Bx=smooth(Bx,Nsm);
By=smooth(By,Nsm);
Bz=smooth(Bz,Nsm);

%options=odeset('AbsTol',10e-1,'RelTol',10e-1,'Stats','off'); 
counter=0
for vxp=-.1:.01:.1
%[t,y]=ode45(@newton,0:.1:2000 ,[24 15 0 vxp .01 .01])
counter=counter+1
load(['autotraject' num2str(counter) '.mat'],'t','y')
%plot(y(:,1),y(:,2))
xgsm=-y(:,1)/Lx*30-15;
zgsm=y(:,2)/Ly*12-9;
subplot(2,2,1)
%plot(y(:,1),y(:,2),y(1,1),y(1,2),'ro',y(end,1),y(end,2),'kx')
hold on
imagesc(xRE,yRE,Bx')
xc=linspace(-45,-15,Nx);
yc=linspace(-9,3,Ny);
[xc,yc]=meshgrid(xc,yc);
hlines=streamslice(xc,yc,-fliplr(Bx')*dy/dx,fliplr(By'));
set(hlines,'Color','w')
set(gca,'xdir','reverse')
plot(xgsm,zgsm,'m',xgsm(1),zgsm(1),'ro',xgsm(end),zgsm(end),'kx')
axis tight
xlim([-30 -15])
xlabel('x/RE')
ylabel('z/RE')
subplot(2,2,2)
%plot(y(:,1),y(:,4),y(1,1),y(1,4),'ro',y(end,1),y(end,4),'kx')
plot(xgsm,y(:,4),xgsm(1),y(1,4),'ro',xgsm(end),y(end,4),'kx')
set(gca,'xdir','reverse')
axis tight
xlabel('x/RE')
ylabel('vxgsm/c')
subplot(2,2,3)
plot(y(:,4),y(:,5),y(1,4),y(1,5),'ro',y(end,4),y(end,5),'kx')
xlabel('vx/c')
ylabel('vzgsm/c')
subplot(2,2,4)
plot(y(:,4),y(:,6),y(1,4),y(1,6),'ro',y(end,4),y(end,6),'kx')
xlabel('vxgsm/c')
ylabel('vygsm/c')
print('-dpdf', ['autotraject' num2str(counter) 'OL.pdf']) 
end
