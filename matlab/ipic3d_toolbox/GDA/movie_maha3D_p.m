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
contours=0;

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

B2=Bx.^2+By.^2+Bz.^2;

file=[dir 'rho_0_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
rho=fread(fid,'real*8');
fclose(fid);
rho=reshape(rho,Nx,Ny,Nz);
rho=squeeze(rho(:,:,iz));


file=[dir 'Pe_par_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Ppar=fread(fid,'real*8');
fclose(fid);
Ppar=reshape(Ppar,Nx,Ny,Nz);
Ppar=squeeze(Ppar(:,:,iz));

file=[dir 'Pe_per1_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Pper1=fread(fid,'real*8');
fclose(fid);
Pper1=reshape(Pper1,Nx,Ny,Nz);
Pper1=squeeze(Pper1(:,:,iz));

file=[dir 'Pe_per2_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Pper2=fread(fid,'real*8');
fclose(fid);
Pper2=reshape(Pper2,Nx,Ny,Nz);
Pper2=squeeze(Pper2(:,:,iz));

Tpar=Ppar./(-rho);
Tper1=Pper1./(-rho);
Tper2=Pper2./(-rho);

adiabatic=Tper1./sqrt(B2);

Ay=zeros(size(Bx));
if (contours)
Ay=vecpot_uniform(xc,yc,Bx*dy/dx,By);
else
xc=linspace(-45,-15,Nx)+dx/2;
yc=linspace(-9,3,Ny)+dy/2;
[xc,yc]=meshgrid(xc,yc);
end

immagine(x,y,Tpar*code_T,['Tepar' ncycle1],[7 18],2,ncycle1, Ygsm)
hold on
if(contours)
%contour(xc,yc,fliplr(Ay'),50,'w')
contour(-15-xc/Lx*30,-9+yc/Ly*12,Ay',50,'w')
%ylim([-5 -1])
%xlim([-35 -20])
else
[sx,sy] = meshgrid(0:Lx/20:Lx,0:Ly/20:Ly);
hlines=streamslice(xc,yc,-fliplr(Bx')*dy/dx,fliplr(By'),2);
set(hlines,'Color','k','linewidth', 1)
set(gcf, 'Renderer', 'zbuffer');
end
name=['Tepar_combo' ncycle1];
print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])
%coplot_uniform(xc,yc,By'*code_B,Ay','x/d_i','y/d_i','By')
%caxis([-5 5])

immagine(x,y,Tper1*code_T,['Teper1' ncycle1],[5 20],0,ncycle1, Ygsm)
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
name=['Teper1_combo' ncycle1];
print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])
%coplot_uniform(xc,yc,By'*code_B,Ay','x/d_i','y/d_i','By')
%caxis([-5 5])


immagine(x,y,Tper2*code_T,['Teper2' ncycle1],[5 20],5,ncycle1, Ygsm)
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
name=['Teper2_combo' ncycle1];
print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])
%coplot_uniform(xc,yc,By'*code_B,Ay','x/d_i','y/d_i','By')
%caxis([-5 5])

immagine(x,y,(Tper1+Tper2)./Tpar./2,['AnisT' ncycle1],[0.5 1.5],2,ncycle1, Ygsm)
immagine(x,y,Tper1./Tper2,['AgyrT' ncycle1],[0.9 1.1],2,ncycle1, Ygsm)
        
immagine(x,y,log10(adiabatic),['Adiabatic' ncycle1],[0 0],5,ncycle1, Ygsm)

end
