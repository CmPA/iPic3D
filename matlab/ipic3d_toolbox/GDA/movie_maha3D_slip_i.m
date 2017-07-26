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


for cycle=000:1000:Ncyc_max

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

Vexbx = (Ey.* Bz - Ez.* By)./B2;
Vexby = (Ez.* Bx - Ex.* Bz)./B2;
Vexbz = (Ex.* By - Ey.* Bx)./B2;

file=[dir 'rho_1_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
rho=fread(fid,'real*8');
fclose(fid);
rho=reshape(rho,Nx,Ny,Nz);
rho=squeeze(rho(:,:,iz));

close all
file=[dir 'Ji_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Vx=fread(fid,'real*8');
fclose(fid);
Vx=reshape(Vx,Nx,Ny,Nz);
Vx=squeeze(Vx(:,:,iz));

file=[dir 'Ji_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Vy=fread(fid,'real*8');
fclose(fid);
Vy=reshape(Vy,Nx,Ny,Nz);
Vy=squeeze(Vy(:,:,iz));

file=[dir 'Ji_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Vz=fread(fid,'real*8');
fclose(fid);
Vz=reshape(Vz,Nx,Ny,Nz);
Vz=squeeze(Vz(:,:,iz));


JdotE = Vx.*Ex+ Vy.*Ey + Vz.*Ez;

Vx=Vx./rho;
Vy=Vy./rho;
Vz=Vz./rho;

VparB= (Vx.*Bx + Vy.*By + Vz.*Bz);

Vparx=VparB.*Bx./B2;
Vpary=VparB.*By./B2;
Vparz=VparB.*Bz./B2;

immagine(x,y,log10(abs(JdotE)),['JdotE' ncycle1],[0 0],0,ncycle1, Ygsm)
immagine(x,y,-(Vx - Vexbx)*code_V,['Vslip,x' ncycle1],[-5 5]*1e2,5,ncycle1, Ygsm)
immagine(x,y,(Vy - Vexby)*code_V,['Vslip,y' ncycle1],[-5 5]*1e2,5,ncycle1, Ygsm)
immagine(x,y,(Vz - Vexbz)*code_V,['Vslip,z' ncycle1],[-5 5]*1e2,5,ncycle1, Ygsm)
immagine(x,y,-(Vx - Vparx - Vexbx)*code_V,['Vslip,perpx' ncycle1],[-5 5]*1e2,5,ncycle1, Ygsm)
immagine(x,y,(Vy - Vpary - Vexby)*code_V,['Vslip,perpy' ncycle1],[-5 5]*1e2,5,ncycle1, Ygsm)
immagine(x,y,(Vz - Vparz - Vexbz)*code_V,['Vslip,perpz' ncycle1],[-5 5]*1e2,5,ncycle1, Ygsm)

hold on

end
