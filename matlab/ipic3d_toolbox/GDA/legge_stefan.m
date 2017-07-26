clear all
close all

dir='/shared02/gianni/StefanO12/data2/'
cycle=50000
%dir='/shared02/gianni/drake3/data/'
%cycle=25000
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
xc=linspace(0, Lx, Nx);
yc=linspace(0, Ly, Ny);


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
if (contours)
Ay=vecpot_uniform(xc,yc,Bx,By);
end


file=[dir 'Je_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Jex=fread(fid,'real*8');
fclose(fid);
Jex=reshape(Jex,Nx,Ny);

file=[dir 'Je_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Jez=fread(fid,'real*8');
fclose(fid);
Jez=reshape(Jez,Nx,Ny);

file=[dir 'Ji_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Jix=fread(fid,'real*8');
fclose(fid);
Jix=reshape(Jix,Nx,Ny);

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

Vex=Jex./rhoe;
Vix=Jix./rhoi;
Nsmooth=10

bxmax=max(abs(Bx(:)));
immagine(xc,yc,Bx,['Bx' ncycle],[-bxmax bxmax],Nsmooth)
hold on
contour(xc,yc,fliplr(Ay'),20,'w')
set(gcf, 'Renderer', 'zbuffer');
name=['Bx_combo' ncycle];
print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])
%coplot_uniform(xc,yc,By'*code_B,Ay','x/d_i','y/d_i','By')
%caxis([-5 5])

bymax=max(abs(By(:)));
immagine(xc,yc,By,['By' ncycle],[-bymax bymax],Nsmooth)
hold on
contour(xc,yc,fliplr(Ay'),20,'w')
set(gcf, 'Renderer', 'zbuffer');
name=['By_combo' ncycle];
print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])


bzmax=max((Bz(:)));
bzmin=min((Bz(:)));
immagine(xc,yc,Bz,['Bz' ncycle],[bzmin bzmax],Nsmooth)
hold on
contour(xc,yc,fliplr(Ay'),20,'w')
set(gcf, 'Renderer', 'zbuffer');
name=['Bz_combo' ncycle];
print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])


vmax=max((Vex(:)));
immagine(xc,yc,Vex,['Vex' ncycle],[-vmax vmax],Nsmooth)
hold on
contour(xc,yc,fliplr(Ay'),20,'w')
set(gcf, 'Renderer', 'zbuffer');
name=['Vex_combo' ncycle];
print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])

vmax=max(abs(Jez(:)));
immagine(xc,yc,Jez,['Jez' ncycle],[-vmax vmax],Nsmooth)
hold on
contour(xc,yc,fliplr(Ay'),20,'w')
set(gcf, 'Renderer', 'zbuffer');
name=['Vex_combo' ncycle];
print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])

vmax=max((Vix(:)));
immagine(xc,yc,Vix,['Vix' ncycle],[-vmax vmax],Nsmooth)
hold on
contour(xc,yc,fliplr(Ay'),20,'w')
set(gcf, 'Renderer', 'zbuffer');
name=['Vix_combo' ncycle];
print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])
