maxp_common

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
%Ygsm=gsmz2y(Lz-zc(iz));


for cycle=Ncyc_ini:1000:Ncyc_max

ncycle=num2str(cycle)
ncycle1=num2str(cycle,'%06d')



global blowup contours
blowup=0;
contours=1;


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

file=[dir 'rho_1_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
rho=fread(fid,'real*8');
fclose(fid);
rho=reshape(rho,Nx,Ny,Nz);
rho=squeeze(rho(:,:,iz));

Vix=Vx./rho;
Viy=Vy./rho;
Viz=Vz./rho;

file=[dir 'Je_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Vx=fread(fid,'real*8');
fclose(fid);
Vx=reshape(Vx,Nx,Ny,Nz);
Vx=squeeze(Vx(:,:,iz));

file=[dir 'Je_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Vy=fread(fid,'real*8');
fclose(fid);
Vy=reshape(Vy,Nx,Ny,Nz);
Vy=squeeze(Vy(:,:,iz));

file=[dir 'Je_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Vz=fread(fid,'real*8');
fclose(fid);
Vz=reshape(Vz,Nx,Ny,Nz);
Vz=squeeze(Vz(:,:,iz));

file=[dir 'rho_0_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
rho=fread(fid,'real*8');
fclose(fid);
rho=reshape(rho,Nx,Ny,Nz);
rho=squeeze(rho(:,:,iz));

Vex=Vx./rho;
Vey=Vy./rho;
Vez=Vz./rho;


file=[dir 'Pe_par_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Ppar=fread(fid,'real*8');
fclose(fid);
Ppar=reshape(Ppar,Nx,Ny,Nz);
Ppar=squeeze(Ppar(:,:,iz));

Tepar=Ppar./(-rho);

vthe=sqrt(Tepar*256);

vlim=[-0.2 0.2];
immagine(x,y,-(Vix-Vex)./vthe,['VRELovthe_l' ncycle1],vlim,3,ncycle1, Ygsm)
immagine(x,y,(Viy-Vey)./vthe,['VRELovthe_n' ncycle1],vlim,3,ncycle1, Ygsm)
immagine(x,y,(Viz-Vez)./vthe,['VRELovthe_m' ncycle1],vlim,3,ncycle1, Ygsm)

vlim=[0.1 0.2];
immagine(x,y,(sqrt((Viz-Vez).^2+(Viy-Vey).^2+(Vix-Vex).^2)./vthe),['VRELovthe_mod' ncycle1],vlim,3,ncycle1, Ygsm)

end
