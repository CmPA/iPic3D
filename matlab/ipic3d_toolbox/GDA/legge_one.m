clear all
close all

dir='/shared/gianni/drake2/part2/'
cycle=25000
%dir='/shared02/gianni/drake4/data/'
%cycle=19000
ncycle=num2str(cycle);

global blowup contours
blowup=0;
contours=1;


addpath '/home/gianni/matlab2/matlab-parsek'
filename=[dir 'settings.hdf'];
Lx=hdf5read(filename,'/collective/Lx'); 
Ly=hdf5read(filename,'/collective/Ly');
Nx=hdf5read(filename,'/collective/Nxc'); 
Ny=hdf5read(filename,'/collective/Nyc');
x=[0 Lx]; 
y=[0 Ly];
xn=linspace(0,Lx,Nx);
yn=linspace(0,Ly,Ny);




file=[dir 'Ji_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
V=reshape(V,Nx,Ny);
file=[dir 'rho_1_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
N=fread(fid,'real*8');
fclose(fid);
N=reshape(N,Nx,Ny);
file=[dir 'B_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Vx=fread(fid,'real*8');
fclose(fid);
Vx=reshape(Vx,Nx,Ny);

V=V./N;

vmin=min(V(:));
vmax=max(V(:));
%immagine_dir(x,y,V,['Jx' ncycle],[vmin vmax],0)
immagine_dir(x,y,V,['Vx' ncycle],[vmin 0],0)
xlim([100 400])
ylim([50 150])

name='Vx'
print('-dpng','-r300',[name '.png'])


vmin=min(Vx(:));
vmax=max(Vx(:));
immagine_dir(x,y,Vx,['By' ncycle],[vmin vmax],0)
xlim([100 400])
ylim([50 150])
name='By'
print('-dpng','-r300',[name '.png'])

close all
xi=100:.1:400;
yi=101*ones(size(xi));
byi=interp2(xn,yn,Vx',xi,yi);
vxi=interp2(xn,yn,V',xi,yi);
subplot(2,1,1)
plot(xi,byi)
title('B_y')
subplot(2,1,2)
plot(xi,-vxi)
title('-V_x')
xlabel('x/d_i')
print('-dpng','satfake.png')
