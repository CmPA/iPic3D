clearvars -except Ncyc_max Ncyc_ini dir fraciz Ygsm
close all

addpath('/home/gianni/Documents/matlab/matlab3')
addpath '/home/gianni/matlab2/matlab-parsek'


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

B0=0.0026

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
y=[-8.7 3.3];

%iz=round(Nz*fraciz);
iz = Nz-round(Nz*(9-Ygsm)/12);
Ygsm=gsmz2y(Lz-zc(iz));



for cycle=Ncyc_ini:1000:Ncyc_max

ncycle=num2str(cycle)
ncycle1=num2str(cycle,'%06d')


global blowup contours
blowup=0;
contours=1;

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

ile=[dir 'E_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Ex=fread(fid,'real*8');
fclose(fid);
Ex=reshape(Ex,Nx,Ny,Nz);
Ex=squeeze(Ex(:,:,iz));

file=[dir 'E_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Ey=fread(fid,'real*8');
fclose(fid);
Ey=reshape(Ey,Nx,Ny,Nz);
Ey=squeeze(Ey(:,:,iz));

file=[dir 'E_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Ez=fread(fid,'real*8');
fclose(fid);
Ez=reshape(Ez,Nx,Ny,Nz);
Ez=squeeze(Ez(:,:,iz));


Ay=zeros(size(Bx));
if (contours)
Ay=vecpot_uniform(xc,yc,Bx*dy/dx,By);
else
xc=linspace(-45,-15,Nx);
yc=linspace(-9,3,Ny);
[xc,yc]=meshgrid(xc,yc);
end


immagine(x,y,-Bx*code_B,['Bl' ncycle1],[-25 25],0,ncycle1, Ygsm)
hold on
if(contours)
%contour(xc,yc,fliplr(Ay'),50,'w')
contour(-15-xc/Lx*30,-9+yc/Ly*12,Ay',50,'w')
%ylim([-5 -1])
%xlim([-35 -20])
else
hlines=streamslice(xc,yc,-fliplr(Bx')*dy/dx,fliplr(By'));
set(hlines,'Color','w')
set(gcf, 'Renderer', 'zbuffer');
end
name=['Bl_combo' ncycle1];
print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])
%coplot_uniform(xc,yc,By'*code_B,Ay','x/d_i','y/d_i','By')
%caxis([-5 5])


immagine(x,y,By*code_B,['Bn' ncycle1],[-10 10],0,ncycle1, Ygsm)
hold on
if(contours)
%contour(xc,yc,fliplr(Ay'),50,'w')
contour(-15-xc/Lx*30,-9+yc/Ly*12,Ay',50,'w')
%ylim([-5 -1])
%xlim([-35 -20])
else
hlines=streamslice(xc,yc,-fliplr(Bx')*dy/dx,fliplr(By'));
set(hlines,'Color','w')
set(gcf, 'Renderer', 'zbuffer');
end
name=['Bn_combo' ncycle1];
print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])
%coplot_uniform(xc,yc,By'*code_B,Ay','x/d_i','y/d_i','By')
%caxis([-5 5])

immagine(x,y,Bz*code_B,['Bm' ncycle1],[-10 10],0,ncycle1, Ygsm)
hold on
if(contours)
%contour(xc,yc,fliplr(Ay'),50,'w')
contour(-15-xc/Lx*30,-9+yc/Ly*12,Ay',50,'w')
%ylim([-5 -1])
%xlim([-35 -20])
else
%starty=linspace(0,Ly,30);
%startx=ones(size(starty))*Ly/4;
%streamline(Bx*dx/dy,By,startx,starty)

hlines=streamslice(xc,yc,-fliplr(Bx')*dy/dx,fliplr(By'));
set(hlines,'Color','w')
end

set(gcf, 'Renderer', 'zbuffer');
name=['Bm_combo' ncycle1];
print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])
%coplot_uniform(xc,yc,By'*code_B,Ay','x/d_i','y/d_i','By')
%caxis([-5 5])

immagine(x,y,-Ex/B0^2,['E_l' ncycle1],[0 0],5,ncycle1, Ygsm)
        hold on
if(contours)
%contour(xc,yc,fliplr(Ay'),50,'w')
contour(-15-xc/Lx*30,-9+yc/Ly*12,Ay',50,'w')
%ylim([-5 -1])
%xlim([-35 -20])
else
%starty=linspace(0,Ly,30);
%startx=ones(size(starty))*Ly/4;
%streamline(Bx*dx/dy,By,startx,starty)

hlines=streamslice(xc,yc,-fliplr(Bx')*dy/dx,fliplr(By'));
set(hlines,'Color','w')
end

set(gcf, 'Renderer', 'zbuffer');
name=['El_combo' ncycle1];
print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])


immagine(x,y,Ez/B0^2,['E_m' ncycle1],[-0.4 0.4],5,ncycle1, Ygsm)
if(contours)
%contour(xc,yc,fliplr(Ay'),50,'w')
contour(-15-xc/Lx*30,-9+yc/Ly*12,Ay',50,'w')
%ylim([-5 -1])
%xlim([-35 -20])
else
%starty=linspace(0,Ly,30);
%startx=ones(size(starty))*Ly/4;
%streamline(Bx*dx/dy,By,startx,starty)

hlines=streamslice(xc,yc,-fliplr(Bx')*dy/dx,fliplr(By'));
set(hlines,'Color','w')
end

set(gcf, 'Renderer', 'zbuffer');
name=['Em_combo' ncycle1];
print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])


immagine(x,y,Ey/B0^2,['E_n' ncycle1],[-0.4 0.4],5,ncycle1, Ygsm)
if(contours)
%contour(xc,yc,fliplr(Ay'),50,'w')
contour(-15-xc/Lx*30,-9+yc/Ly*12,Ay',50,'w')
%ylim([-5 -1])
%xlim([-35 -20])
else
%starty=linspace(0,Ly,30);
%startx=ones(size(starty))*Ly/4;
%streamline(Bx*dx/dy,By,startx,starty)

hlines=streamslice(xc,yc,-fliplr(Bx')*dy/dx,fliplr(By'));
set(hlines,'Color','w')
end

set(gcf, 'Renderer', 'zbuffer');
name=['En_combo' ncycle1];
print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])
end
