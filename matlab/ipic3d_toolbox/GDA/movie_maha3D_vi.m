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

Vx=Vx./rho;
Vy=Vy./rho;
Vz=Vz./rho;


vlim=[-500 500];
immagine(x,y,-Vx*code_V,['Vl' ncycle1],vlim,3,ncycle1, Ygsm)
immagine(x,y,Vy*code_V,['Vn' ncycle1],vlim,3,ncycle1, Ygsm)
immagine(x,y,Vz*code_V,['Vm' ncycle1],vlim,3,ncycle1, Ygsm)

end
