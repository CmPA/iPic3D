clear all
close all

for ix=150:250
nix=num2str(ix,'%06d')
%for cycle=0:1000:63000
for cycle=63000:63000
ncycle=num2str(cycle)
ncycle1=num2str(cycle,'%06d')

dir='/data/gianni/gda/HRmaha3D1/'


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
x=[-15 -45]; 
y=[-9 3];
z=[-10 2];

%ix=round(2*Nx/3)
%ix=230

close all
file=[dir 'Je_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Bx=fread(fid,'real*8');
fclose(fid);
Bx=reshape(Bx,Nx,Ny,Nz);
Bx=squeeze(Bx(ix,:,:));

file=[dir 'Je_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
By=fread(fid,'real*8');
fclose(fid);
By=reshape(By,Nx,Ny,Nz);
By=squeeze(By(ix,:,:));

file=[dir 'Je_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Bz=fread(fid,'real*8');
fclose(fid);
Bz=reshape(Bz,Nx,Ny,Nz);
Bz=squeeze(Bz(ix,:,:));

skippa=0
if(skippa)

Ay=zeros(size(Bx));
if (contours)
Ay=vecpot_uniform(xc,yc,Bx*dy/dx,By);
else
xc=linspace(-45,-15,Nx);
yc=linspace(-9,3,Ny);
[xc,yc]=meshgrid(xc,yc);
end


immagine(z,y,Bx'*code_J,['Jl' ncycle1],[-4 1],5,ncycle1,1,'z/R_E','y/R_E')
%colormap(autumn)
ylim([-3.5 -0])
hold on
if(contours)
%contour(xc,yc,fliplr(Ay'),50,'w')
%contour(-15-xc/Lx*30,-9+yc/Ly*12,Ay',50,'w')
%ylim([-5 -1])
%xlim([-35 -20])
else
hlines=streamslice(xc,yc,-fliplr(Bx')*dy/dx,fliplr(By'));
set(hlines,'Color','w')
set(gcf, 'Renderer', 'zbuffer');
end
name=['Jl_combo' ncycle1];
print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])
%coplot_uniform(xc,yc,By'*code_B,Ay','x/d_i','y/d_i','By')
%caxis([-5 5])


immagine(z,y,By'*code_J,['Jn' ncycle1],[-1 1],5,ncycle1,2,'z/R_E','y/R_E')
%colormap(autumn)
ylim([-3.5 -0])
hold on
if(contours)
%contour(xc,yc,fliplr(Ay'),50,'w')
%contour(-15-xc/Lx*30,-9+yc/Ly*12,Ay',50,'w')
%ylim([-5 -1])
%xlim([-35 -20])
else
hlines=streamslice(xc,yc,-fliplr(Bx')*dy/dx,fliplr(By'));
set(hlines,'Color','w')
set(gcf, 'Renderer', 'zbuffer');
end
name=['Jn_combo' ncycle1];
print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])
%coplot_uniform(xc,yc,By'*code_B,Ay','x/d_i','y/d_i','By')
%caxis([-5 5])

end

immagine_yz(z,y,Bz'*code_J,['Jm' nix],[-1 4],5,ncycle1, ix*dx,3,'z/R_E','y/R_E')
%colormap(autumn)
ylim([-3.5 -0])
hold on
if(contours)
%contour(xc,yc,fliplr(Ay'),50,'w')
%contour(-15-xc/Lx*30,-9+yc/Ly*12,Ay',50,'w')
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
name=['Jm_combo' nix];
print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])
%coplot_uniform(xc,yc,By'*code_B,Ay','x/d_i','y/d_i','By')
%caxis([-5 5])
end
end
