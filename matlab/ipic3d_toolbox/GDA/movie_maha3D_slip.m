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
iz=iz-10:iz+10;
%Ygsm=gsmz2y(Lz-zc(iz),3));


for cycle=Ncyc_ini:1000:Ncyc_max

ncycle=num2str(cycle)
ncycle1=num2str(cycle,'%06d')



global blowup contours
blowup=0;
contours=1;

close all
file=[dir 'E_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Ex=fread(fid,'real*8');
fclose(fid);
Ex=reshape(Ex,Nx,Ny,Nz);

file=[dir 'E_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Ey=fread(fid,'real*8');
fclose(fid);
Ey=reshape(Ey,Nx,Ny,Nz);


file=[dir 'E_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Ez=fread(fid,'real*8');
fclose(fid);
Ez=reshape(Ez,Nx,Ny,Nz);

Ex=squeeze(mean(Ex(:,:,iz),3));
Ey=squeeze(mean(Ey(:,:,iz),3));
Ez=squeeze(mean(Ez(:,:,iz),3));

file=[dir 'B_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Bx=fread(fid,'real*8');
fclose(fid);
Bx=reshape(Bx,Nx,Ny,Nz);


file=[dir 'B_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
By=fread(fid,'real*8');
fclose(fid);
By=reshape(By,Nx,Ny,Nz);


file=[dir 'B_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Bz=fread(fid,'real*8');
fclose(fid);
Bz=reshape(Bz,Nx,Ny,Nz);

[X Y Z] = meshgrid(dx/2:dx:Lx-dx/2,dy/2:dy:Ly-dy/2,dz/2:dz:Lz-dz/2);
size(permute(Bx,[2 1 3]))
size(X)

[Jx, Jy, Jz]=curl(X,Y,Z,permute(Bx,[2 1 3]),permute(By,[2 1 3]),permute(Bz,[2 1 3])); % Valid only because dx=dy=dz
Jx=permute(Jx,[2 1 3])/4/pi;
Jy=permute(Jy,[2 1 3])/4/pi;
Jz=permute(Jz,[2 1 3])/4/pi;


Bx=squeeze(mean(Bx(:,:,iz),3));
By=squeeze(mean(By(:,:,iz),3));
Bz=squeeze(mean(Bz(:,:,iz),3));

Jx=squeeze(mean(Jx(:,:,iz),3));
Jy=squeeze(mean(Jy(:,:,iz),3));
Jz=squeeze(mean(Jz(:,:,iz),3));

B2=Bx.^2+By.^2+Bz.^2;

Vexbx = (Ey.* Bz - Ez.* By)./B2;
Vexby = (Ez.* Bx - Ex.* Bz)./B2;
Vexbz = (Ex.* By - Ey.* Bx)./B2;

file=[dir 'rho_0_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
rho=fread(fid,'real*8');
fclose(fid);
rho=reshape(rho,Nx,Ny,Nz);
rho=squeeze(mean(rho(:,:,iz),3));

file=[dir 'rho_1_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
rhoi=fread(fid,'real*8');
fclose(fid);
rhoi=reshape(rhoi,Nx,Ny,Nz);
rhoi=squeeze(mean(rhoi(:,:,iz),3));


rhoc=rhoi+rho;

close all
file=[dir 'Je_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Vx=fread(fid,'real*8');
fclose(fid);
Vx=reshape(Vx,Nx,Ny,Nz);
Vx=squeeze(mean(Vx(:,:,iz),3));

file=[dir 'Je_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Vy=fread(fid,'real*8');
fclose(fid);
Vy=reshape(Vy,Nx,Ny,Nz);
Vy=squeeze(mean(Vy(:,:,iz),3));

file=[dir 'Je_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Vz=fread(fid,'real*8');
fclose(fid);
Vz=reshape(Vz,Nx,Ny,Nz);
Vz=squeeze(mean(Vz(:,:,iz),3));


JedotE = Vx.*Ex+ Vy.*Ey + Vz.*Ez;

JdotE = Jx.*Ex+ Jy.*Ey + Jz.*Ez;


Vx=Vx./rho;
Vy=Vy./rho;
Vz=Vz./rho;

Zx = Ex + (Vy.* Bz - Vz.* By);
Zy = Ey + (Vz.* Bx - Vx.* Bz);
Zz = Ez + (Vx.* By - Vy.* Bx);
Zenitani = (Jx.*Zx+Jy.*Zy+Jz.*Zz) -rhoc.*(Vx.*Ex+Vy.*Ey+Vz.*Ez);


VparB= (Vx.*Bx + Vy.*By + Vz.*Bz);

Vparx=VparB.*Bx./B2;
Vpary=VparB.*By./B2;
Vparz=VparB.*Bz./B2;


immagine(x,y,JdotE,['JdotE' ncycle1],[-2 2]*1e-10,6,ncycle1, Ygsm)
immagine(x,y,Zenitani,['Zenitani' ncycle1],[-1 1]*1e-10,6,ncycle1, Ygsm)


immagine(x,y,-(Vx)*code_V,['Ve,x' ncycle1],[-1 1]*1e3,3,ncycle1, Ygsm)
immagine(x,y,(Vy)*code_V,['Ve,y' ncycle1],[-1 1]*1e3,3,ncycle1, Ygsm)
immagine(x,y,(Vz)*code_V,['Ve,z' ncycle1],[-1 1]*1e3,3,ncycle1, Ygsm)



Vslip=sqrt((Vx - Vexbx).^2+(Vy - Vexby).^2+ (Vz - Vexbz).^2);

immagine(x,y,-(Vx - Vexbx)*code_V,['Vslip,x' ncycle1],[-2 2]*1e3,5,ncycle1, Ygsm)
immagine(x,y,(Vy - Vexby)*code_V,['Vslip,y' ncycle1],[-2 2]*1e3,5,ncycle1, Ygsm)
immagine(x,y,(Vz - Vexbz)*code_V,['Vslip,z' ncycle1],[-2 2]*1e3,5,ncycle1, Ygsm)

immagine(x,y,Vslip*code_V,['Vslip' ncycle1],[-2 2]*1e3,5,ncycle1, Ygsm)

Vslipperp = sqrt( (Vx - Vparx - Vexbx).^2 + (Vy - Vpary - Vexby).^2 + (Vz - Vparz - Vexbz).^2  );

immagine(x,y,-(Vx - Vparx - Vexbx)*code_V,['Vslip,perpx' ncycle1],[-2 2]*1e3,5,ncycle1, Ygsm)
immagine(x,y,(Vy - Vpary - Vexby)*code_V,['Vslip,perpy' ncycle1],[-2 2]*1e3,5,ncycle1, Ygsm)
immagine(x,y,(Vz - Vparz - Vexbz)*code_V,['Vslip,perpz' ncycle1],[-2 2]*1e3,5,ncycle1, Ygsm)
immagine(x,y,Vslipperp*code_V,['Vslip,perp' ncycle1],[-2 2]*1e3,5,ncycle1, Ygsm)

hold on

end
