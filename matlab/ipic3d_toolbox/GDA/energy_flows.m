maxp_common


iz = Nz-round(Nz*(max(Zgsmrange)-Ygsm)/(max(Zgsmrange)-min(Zgsmrange)));
Ygsm=gsmz2y(Lz-zc(iz));



for cycle=Ncyc_ini:1000:Ncyc_max

ncycle=num2str(cycle)
ncycle1=num2str(cycle,'%06d')


global color_choice symmetric_color
blowup=0;
contours=1;

close all
file=[dir 'B_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
Bx=reshape(V,Nx,Ny,Nz);
BxXY=squeeze(Bx(:,:,iz));

file=[dir 'B_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
By=reshape(V,Nx,Ny,Nz);
ByXY=squeeze(By(:,:,iz));

file=[dir 'B_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
Bz=reshape(V,Nx,Ny,Nz);
BzXY=squeeze(Bz(:,:,iz));


ile=[dir 'E_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
Ex=reshape(V,Nx,Ny,Nz);
ExXY=squeeze(Ex(:,:,iz));


file=[dir 'E_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
Ey=reshape(V,Nx,Ny,Nz);
EyXY=squeeze(Ey(:,:,iz));

file=[dir 'E_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
Ez=reshape(V,Nx,Ny,Nz);
EzXY=squeeze(Ez(:,:,iz));

file=[dir 'Ji_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
Jix=reshape(V,Nx,Ny,Nz);
JixXY=squeeze(Jix(:,:,iz));

file=[dir 'Ji_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
Jiy=reshape(V,Nx,Ny,Nz);
JiyXY=squeeze(Jiy(:,:,iz));


file=[dir 'Ji_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
Jiz=reshape(V,Nx,Ny,Nz);
JizXY=squeeze(Jiz(:,:,iz));


file=[dir 'Je_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
Jex=reshape(V,Nx,Ny,Nz);
JexXY=squeeze(Jex(:,:,iz));

file=[dir 'Je_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
Jey=reshape(V,Nx,Ny,Nz);
JeyXY=squeeze(Jey(:,:,iz));


file=[dir 'Je_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
Jez=reshape(V,Nx,Ny,Nz);
JezXY=squeeze(Jez(:,:,iz));


[SxXY,SyXY,SzXY]=cross_prod(ExXY,EyXY,EzXY,BxXY,ByXY,BzXY);
SxXY=SxXY/4/pi;
SyXY=SyXY/4/pi;
SzXY=SzXY/4/pi;

[Sx,Sy,Sz]=cross_prod(Ex,Ey,Ez,Bx,By,Bz);
SxXZ=sum(Sx,2)*dy/4/pi;
SyXZ=sum(Sy,2)*dy/4/pi;
SzXZ=sum(Sz,2)*dy/4/pi;

[X Y Z] = meshgrid(dx/2:dx:Lx-dx/2,dy/2:dy:Ly-dy/2,dz/2:dz:Lz-dz/2);

divS=divergence(X,Y,Z,permute(Sx,[2 1 3]),permute(Sy,[2 1 3]),permute(Sz,[2 1 3])); % Valid only because dx=dy=dz
divS=permute(divS,[2 1 3])/4/pi;

JidotEXY=JixXY.*ExXY+JiyXY.*EyXY+JizXY.*EzXY;
JidotEXZ=sum(Jix.*Ex+Jiy.*Ey+Jiz.*Ez,2)*dy;

JedotEXY=JexXY.*ExXY+JeyXY.*EyXY+JezXY.*EzXY;
JedotEXZ=sum(Jex.*Ex+Jey.*Ey+Jez.*Ez,2)*dy;

divSXY=squeeze(divS(:,:,iz));
divSXZ=sum(divS,2)*dy;

Ay=zeros(size(Bx));
if (contours)
Ay=vecpot_uniform(xc,yc,BxXY*dy/dx,ByXY);
else
xc=linspace(-45,-15,Nx);
yc=linspace(-3,9,Ny);
[xc,yc]=meshgrid(xc,yc);
end

global color_choice symmetric_color
symmetric_color=1 %1 yes (max value) -1 yes (min value) 0 no
color_choice=1 % 1 is David's color scale, 2 is hot 0 is jet
%
%	Sl
%
close all
symmetric_color=-1
ncontours=30
contour_color= 'g'
mu0=1.25663706e-6 ;
code_S=code_B*code_E*4*pi/mu0;
immagine(fliplr(Xgsmrange),Zgsmrange,-SxXY,['SlXY' ncycle1],[-1 1]*0e5,0,ncycle1, Ygsm)
hold on
if(contours)
contour(gsmx(xc),gsmy2z(yc),Ay',ncontours,'color',contour_color);
else
hlines=streamslice(xc,yc,-fliplr(BxXY')*dy/dx,fliplr(ByXY'));
set(hlines,'Color','g')
set(gcf, 'Renderer', 'zbuffer');
end
name=['Sl_combo' ncycle1];
%print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),-SxXZ,['SlXZ' ncycle1],[-1 1]*0e6,0,ncycle1,[],3,'x/R_E','y/R_E','Sl [1E-12] W/m^2');

others=0
if(others)
%
%	Sm
%
close all
symmetric_color=-1 
mu0=1.25663706e-6 ;
code_S=code_B*code_E*4*pi/mu0;
immagine(fliplr(Xgsmrange),Zgsmrange,SzXY,['SmXY' ncycle1],[-1 1]*0e5,0,ncycle1, Ygsm)
hold on
if(contours)
contour(gsmx(xc),gsmy2z(yc),Ay',ncontours,'color',contour_color)
else
hlines=streamslice(xc,yc,-fliplr(BxXY')*dy/dx,fliplr(ByXY'));
set(hlines,'Color','w')
set(gcf, 'Renderer', 'zbuffer');
end
name=['Sm_combo' ncycle1];
%print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),SzXZ,['SmXZ' ncycle1],[-3 3]*0e7,0,ncycle1,[],3,'x/R_E','y/R_E','Sm');

%
%	Sn
%
close all
symmetric_color=1 
mu0=1.25663706e-6 ;
code_S=code_B*code_E*4*pi/mu0;
immagine(fliplr(Xgsmrange),Zgsmrange,SyXY,['SnXY' ncycle1],[-1 1]*0e5,0,ncycle1, Ygsm)
hold on
if(contours)
contour(gsmx(xc),gsmy2z(yc),Ay',ncontours,'color',contour_color)
else
hlines=streamslice(xc,yc,-fliplr(BxXY')*dy/dx,fliplr(ByXY'));
set(hlines,'Color','w')
set(gcf, 'Renderer', 'zbuffer');
end
name=['Sn_combo' ncycle1];
%print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),SyXZ,['SnXZ' ncycle1],[-3 3]*0e7,0,ncycle1,[],3,'x/R_E','y/R_E','Sn');

end

%
%	Ji.E
%
close all
symmetric_color=1 
immagine(fliplr(Xgsmrange),Zgsmrange,JidotEXY,['JidotEXY' ncycle1],[-1 1]*0e5,3,ncycle1, Ygsm)
hold on
if(contours)
contour(gsmx(xc),gsmy2z(yc),Ay',ncontours,'color',contour_color)
else
hlines=streamslice(xc,yc,-fliplr(BxXY')*dy/dx,fliplr(ByXY'));
set(hlines,'Color','w')
set(gcf, 'Renderer', 'zbuffer');
end
name=['JidotE_combo' ncycle1];
%print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),JidotEXZ,['JidotEXZ' ncycle1],[-1 1]*2e-6,3,ncycle1,[],3,'x/R_E','y/R_E','JidotE [1E-12] W/m^2');

%
%	Je.E
%
close all
symmetric_color=1 
immagine(fliplr(Xgsmrange),Zgsmrange,JedotEXY,['JedotEXY' ncycle1],[-1 1]*0e5,3,ncycle1, Ygsm)
hold on
if(contours)
contour(gsmx(xc),gsmy2z(yc),Ay',ncontours,'color',contour_color)
else
hlines=streamslice(xc,yc,-fliplr(BxXY')*dy/dx,fliplr(ByXY'));
set(hlines,'Color','w')
set(gcf, 'Renderer', 'zbuffer');
end
name=['JedotE_combo' ncycle1];
%print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),JedotEXZ,['JedotEXZ' ncycle1],[-1 1]*2e-6,3,ncycle1,[],3,'x/R_E','y/R_E','JedotE [1E-12] W/m^2');



%
%	J.E
%
close all
symmetric_color=1
color_choice=3 
immagine(fliplr(Xgsmrange),Zgsmrange,JedotEXY+JidotEXY,['JdotEXY' ncycle1],[-1 1]*0e5,3,ncycle1, Ygsm)
hold on
if(contours)
contour(gsmx(xc),gsmy2z(yc),Ay',ncontours,'color',contour_color)
else
hlines=streamslice(xc,yc,-fliplr(BxXY')*dy/dx,fliplr(ByXY'));
set(hlines,'Color','w')
set(gcf, 'Renderer', 'zbuffer');
end
name=['JdotE_combo' ncycle1];
%print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),JedotEXZ+JidotEXZ,['JdotEXZ' ncycle1],[-1 1]*.5e-6,6,ncycle1,[],3,'x/R_E','y/R_E','JdotE [1E-12] W/m^2');



%
%	divS
%
close all
symmetric_color=-1
color_choice=3 
immagine(fliplr(Xgsmrange),Zgsmrange,divSXY,['divSXY' ncycle1],[-1 1]*0e5,3,ncycle1, Ygsm)
hold on
if(contours)
contour(gsmx(xc),gsmy2z(yc),Ay',ncontours,'color',contour_color)
else
hlines=streamslice(xc,yc,-fliplr(BxXY')*dy/dx,fliplr(ByXY'));
set(hlines,'Color','w')
set(gcf, 'Renderer', 'zbuffer');
end
name=['divS_combo' ncycle1];
%print('-depsc','-r300',[name '.eps'])
print('-dpng','-r300',[name '.png'])

immagine_xy(gsmx([0 Lx]),gsmz2y([0 Lz]),divSXZ,['divSXZ' ncycle1],[-1 1]*.5e-6,3,ncycle1,[],3,'x/R_E','y/R_E','divS [1E-12] W/m^2');





end
