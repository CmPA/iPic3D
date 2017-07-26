close all
addpath(genpath('../../ipic3d_toolbox'));
dir='/Users/gianni/Desktop/tred68/';


iz=50;
NNcyc=25
for cycle=0:1000:NNcyc*1000

leggo=1;
if(leggo==1)

% [Bx3d,By3d,Bz3d,Nx,Ny,Nz]=read_binVTK_vector(dir,'B',cycle);
% Bx=squeeze(Bx3d(:,:,iz));
% By=squeeze(By3d(:,:,iz));
% Bz=squeeze(Bz3d(:,:,iz));
% [Ex3d,Ey3d,Ez3d,Nx,Ny,Nz]=read_binVTK_vector(dir,'E',cycle);
% Ex=squeeze(Ex3d(:,:,iz));
% Ey=squeeze(Ey3d(:,:,iz));
% Ez=squeeze(Ez3d(:,:,iz));
% [Jex3d,Jey3d,Jez3d,Nx,Ny,Nz]=read_binVTK_vector(dir,'Je',cycle);
% Jex=squeeze(Jex3d(:,:,iz));
% Jey=squeeze(Jey3d(:,:,iz));
% Jez=squeeze(Jez3d(:,:,iz));
% [Jix3d,Jiy3d,Jiz3d,Nx,Ny,Nz]=read_binVTK_vector(dir,'Ji',cycle);
% Jix=squeeze(Jix3d(:,:,iz));
% Jiy=squeeze(Jiy3d(:,:,iz));
% Jiz=squeeze(Jiz3d(:,:,iz));
% 
% [Az3d,Nx,Ny,Nz]=read_binVTK_scalar(dir,'Az',cycle);
% Az=squeeze(Az3d(:,:,iz));

[Ag3d,Nx,Ny,Nz,dx,dy,dz]=read_binVTK_scalar(dir,'Agyro_aunai_shift',cycle);
ixplt=round(Nx/2);
iyplt=round(Ny/2);
izplt=round(Nz/2);
Lx=dx*Nx;
Ly=dy*Ny;
Lz=dz*Nz;
xv=linspace(0,Lx,Nx);
yv=linspace(0,Ly,Ny);
zv=linspace(0,Lz,Nz);
Ncut=3;
xv=xv(Ncut:end-Ncut);
yv=yv(Ncut:end-Ncut);
zv=zv(Ncut:end-Ncut);
AgXY=squeeze(Ag3d(Ncut:end-Ncut,Ncut:end-Ncut,izplt));
AgYZ=squeeze(Ag3d(ixplt,Ncut:end-Ncut,Ncut:end-Ncut));
AgXZ=squeeze(Ag3d(Ncut:end-Ncut,iyplt,Ncut:end-Ncut));
% 
% [rhoe3d,rhoi3d,Nx,Ny,Nz]=read_binVTK_multiscalar(dir,'rho',cycle);
% rhoe=squeeze(rhoe3d(:,:,iz));
% rhoi=squeeze(rhoi3d(:,:,iz));
% 
% [Pixx,Pixy,Pixz,Piyy,Piyz,Pizz,Pipar,Piper1,Piper2,Pieps,Nx,Ny,Nz] = read_binVTK_pressure(dir,'Pi',cycle);
% [Pexx,Pexy,Pexz,Peyy,Peyz,Pezz,Pepar,Peper1,Peper2,Peeps,Nx,Ny,Nz] = read_binVTK_pressure(dir,'Pe',cycle);
% 
% B=sqrt(Bx.*Bx+By.*By+Bz.*Bz);
% Epar=(Ex.*Bx+Ey.*By+Ez.*Bz)./B;

% agyro
end

figure(1)
imagesc(xv,yv,AgXY')
xlabel('x/d_i','fontsize',[15])
ylabel('y/d_i','fontsize',[15])
set(gca,'fontsize',[15])
figure(2)
imagesc(xv,zv,AgXZ')
xlabel('x/d_i','fontsize',[15])
ylabel('z/d_i','fontsize',[15])
set(gca,'fontsize',[15])
%figure(3)
% imagesc(yv,zv,AgYZ')
% xlabel('y')
% ylabel('z')
pause
end
