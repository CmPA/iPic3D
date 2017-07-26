maxp_common

%xc=linspace(-45, -15, Nx);
%yc=linspace(-9, 3, Ny);
xc=linspace(0, Lx, Nx);
yc=linspace(0, Ly, Ny);
zc=linspace(0, Lz, Nz);
x=[-15 -45];
y=[-8.7 3.3];

%iz=round(Nz*fraciz);
iz = Nz-round(Nz*(9-Ygsm)/12);
%Ygsm=gsmz2y(Lz-zc(iz));


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

file=[dir 'rho_1_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
rho=fread(fid,'real*8');
fclose(fid);
rho=reshape(rho,Nx,Ny,Nz);
rho=squeeze(rho(:,:,iz));


file=[dir 'Pi_xx_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Ppar=fread(fid,'real*8');
fclose(fid);
Ppar=reshape(Ppar,Nx,Ny,Nz);
Ppar=squeeze(Ppar(:,:,iz));

file=[dir 'Pi_per1_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Pper1=fread(fid,'real*8');
fclose(fid);
Pper1=reshape(Pper1,Nx,Ny,Nz);
Pper1=squeeze(Pper1(:,:,iz));

file=[dir 'Pi_per2_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Pper2=fread(fid,'real*8');
fclose(fid);
Pper2=reshape(Pper2,Nx,Ny,Nz);
Pper2=squeeze(Pper2(:,:,iz));

Tper=(Pper1+Pper2)./(rho);

file=[dir 'Pi_yy_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Pper1=fread(fid,'real*8');
fclose(fid);
Pper1=reshape(Pper1,Nx,Ny,Nz);
Pper1=squeeze(Pper1(:,:,iz));

file=[dir 'Pi_zz_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Pper2=fread(fid,'real*8');
fclose(fid);
Pper2=reshape(Pper2,Nx,Ny,Nz);
Pper2=squeeze(Pper2(:,:,iz));

TXX=Ppar./(rho);
TYY=Pper1./(rho);
TZZ=Pper2./(rho);

work=0
if(work)
file=[dir 'WorkPi_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
WorkPi_x=fread(fid,'real*8');
fclose(fid);
WorkPi_x=reshape(WorkPi_x,Nx,Ny,Nz);
WorkPi_x=squeeze(WorkPi_x(:,:,iz));

file=[dir 'WorkPi_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
WorkPi_y=fread(fid,'real*8');
fclose(fid);
WorkPi_y=reshape(WorkPi_y,Nx,Ny,Nz);
WorkPi_y=squeeze(WorkPi_y(:,:,iz));

file=[dir 'WorkPi_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
WorkPi_z=fread(fid,'real*8');
fclose(fid);
WorkPi_z=reshape(WorkPi_z,Nx,Ny,Nz);
WorkPi_z=squeeze(WorkPi_z(:,:,iz));

file=[dir 'Enthi_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Enthi_x=fread(fid,'real*8');
fclose(fid);
Enthi_x=reshape(Enthi_x,Nx,Ny,Nz);
Enthi_x=squeeze(Enthi_x(:,:,iz));

file=[dir 'Enthi_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Enthi_y=fread(fid,'real*8');
fclose(fid);
Enthi_y=reshape(Enthi_y,Nx,Ny,Nz);
Enthi_y=squeeze(Enthi_y(:,:,iz));

file=[dir 'Enthi_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Enthi_z=fread(fid,'real*8');
fclose(fid);
Enthi_z=reshape(Enthi_z,Nx,Ny,Nz);
Enthi_z=squeeze(Enthi_z(:,:,iz));

file=[dir 'PxBi_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
PxBi_x=fread(fid,'real*8');
fclose(fid);
PxBi_x=reshape(PxBi_x,Nx,Ny,Nz);
PxBi_x=squeeze(PxBi_x(:,:,iz));

file=[dir 'PxBi_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
PxBi_y=fread(fid,'real*8');
fclose(fid);
PxBi_y=reshape(PxBi_y,Nx,Ny,Nz);
PxBi_y=squeeze(PxBi_y(:,:,iz));

file=[dir 'PxBi_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
PxBi_z=fread(fid,'real*8');
fclose(fid);
PxBi_z=reshape(PxBi_z,Nx,Ny,Nz);
PxBi_z=squeeze(PxBi_z(:,:,iz));

end

Ay=zeros(size(Bx));
if (contours)
Ay=vecpot_uniform(xc,yc,Bx*dy/dx,By);
else
xc=linspace(-45,-15,Nx);
yc=linspace(-9,3,Ny);
[xc,yc]=meshgrid(xc,yc);
end

adiabatic=(Tper)./sqrt(Bx.^2+By.^2+Bz.^2);

immagine(x,y,TXX*code_T,['XZTiX' ncycle1],[30 80],3,ncycle1, Ygsm)
hold on
if(contours)
%contour(xc,yc,fliplr(Ay'),50,'w')
contour(-15-xc/Lx*30,-9+yc/Ly*12,Ay',20,'w')
%ylim([-5 -1])
%xlim([-35 -20])
else
hlines=streamslice(xc,yc,-fliplr(Bx')*dy/dx,fliplr(By'));
set(hlines,'Color','w')
set(gcf, 'Renderer', 'zbuffer');
end
name=['Tipar_combo' ncycle1];
print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])
%coplot_uniform(xc,yc,By'*code_B,Ay','x/d_i','y/d_i','By')
%caxis([-5 5])

immagine(x,y,TYY*code_T,['XZTiY' ncycle1],[30 70],3,ncycle1, Ygsm)
hold on
if(contours)
%contour(xc,yc,fliplr(Ay'),50,'w')
contour(-15-xc/Lx*30,-9+yc/Ly*12,Ay',20,'w')
%ylim([-5 -1])
%xlim([-35 -20])
else
hlines=streamslice(xc,yc,-fliplr(Bx')*dy/dx,fliplr(By'));
set(hlines,'Color','w')
set(gcf, 'Renderer', 'zbuffer');
end
name=['Tiper1_combo' ncycle1];
print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])
%coplot_uniform(xc,yc,By'*code_B,Ay','x/d_i','y/d_i','By')
%caxis([-5 5])


immagine(x,y,TZZ*code_T,['XZTiZ' ncycle1],[30 70],3,ncycle1, Ygsm)
hold on
if(contours)
%contour(xc,yc,fliplr(Ay'),50,'w')
contour(-15-xc/Lx*30,-9+yc/Ly*12,Ay',20,'w')
%ylim([-5 -1])
%xlim([-35 -20])
else
hlines=streamslice(xc,yc,-fliplr(Bx')*dy/dx,fliplr(By'));
set(hlines,'Color','w')
set(gcf, 'Renderer', 'zbuffer');
end
name=['Tiper2_combo' ncycle1];
print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])
%coplot_uniform(xc,yc,By'*code_B,Ay','x/d_i','y/d_i','By')
%caxis([-5 5])

if(work)
nsmooth=8;
val_range=[-.5 .5]*1e-9;
immagine(x,y,-Enthi_x,['XZEnthi_x' ncycle1],val_range,nsmooth,ncycle1, Ygsm)        
immagine(x,y,-Enthi_y,['XZEnthi_y' ncycle1],val_range,nsmooth,ncycle1, Ygsm) 
immagine(x,y,-Enthi_z,['XZEnthi_z' ncycle1],val_range,nsmooth,ncycle1, Ygsm) 

immagine(x,y,WorkPi_x,['XZWorkPi_x' ncycle1],val_range,nsmooth,ncycle1, Ygsm)  
immagine(x,y,WorkPi_y,['XZWorkPi_y' ncycle1],val_range,nsmooth,ncycle1, Ygsm)
immagine(x,y,WorkPi_z,['XZWorkPi_z' ncycle1],val_range,nsmooth,ncycle1, Ygsm) 
   
immagine(x,y,PxBi_x,['XZPxBi_x' ncycle1],val_range,nsmooth,ncycle1, Ygsm) 
immagine(x,y,PxBi_y,['XZPxBi_y' ncycle1],val_range,nsmooth,ncycle1, Ygsm) 
immagine(x,y,PxBi_z,['XZPxBi_z' ncycle1],val_range,nsmooth,ncycle1, Ygsm) 

val_range=[0.0 0.0 ];
immagine(x,y,log10(adiabatic),['XZAdiab' ncycle1],val_range,nsmooth,ncycle1, Ygsm) 
end
        
end
