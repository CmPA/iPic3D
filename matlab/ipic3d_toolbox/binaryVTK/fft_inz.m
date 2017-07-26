close all
addpath(genpath('../../ipic3d_toolbox'));
dir='/shared/gianni/tred70/';

for cycle=21000:1000:21000

leggo=1;
if(leggo==1)

[Bx,By,Bz,Nx,Ny,Nz]=read_binVTK_vector(dir,'B',cycle);
[Ex,Ey,Ez,Nx,Ny,Nz]=read_binVTK_vector(dir,'E',cycle);
% [Jex,Jey,Jez,Nx,Ny,Nz]=read_binVTK_vector(dir,'Je',cycle);
% [Jix,Jiy,Jiz,Nx,Ny,Nz]=read_binVTK_vector(dir,'Ji',cycle);
% 
% [Az,Nx,Ny,Nz]=read_binVTK_scalar(dir,'Az',cycle);
% [rhoe,rhoi,Nx,Ny,Nz]=read_binVTK_multiscalar(dir,'rho',cycle);
% [Pixx,Pixy,Pixz,Piyy,Piyz,Pizz,Pipar,Piper1,Piper2,Pieps,Nx,Ny,Nz] = read_binVTK_pressure(dir,'Pi',cycle);
% [Pexx,Pexy,Pexz,Peyy,Peyz,Pezz,Pepar,Peper1,Peper2,Peeps,Nx,Ny,Nz] = read_binVTK_pressure(dir,'Pe',cycle);
% 
 B=sqrt(Bx.*Bx+By.*By+Bz.*Bz);
 Epar=(Ex.*Bx+Ey.*By+Ez.*Bz)./B;
% 
% Te=(Pexx+Peyy+Pezz)./(-rhoe);
% Ti=(Pixx+Piyy+Pizz)./rhoi;
end

Lx=40;Ly=20;Lz=10;

dx=Lx/Nx;
dy=Ly/Ny;
dz=Lz/Nz;

[X Y] = meshgrid(dx/2:dx:Lx-dx/2,dy/2:dy:Ly-dy/2);

[X Y] = meshgrid(0:dx:Lx-dx,0:dy:Ly-dy);

spectrum = fft(Epar(:,70:end-70,:),[],3);
h=figure(1)
for i=1:9
subplot(3,3,i)
imagesc(real(spectrum(:,:,i))')
xlabel('x/d_i')
ylabel('y/d_i')
title(['E_{||}(x,y) m_z=' num2str(i-1) ])
caxis([-2 2]*1e-3)
end
saveas(h,'tred70.fig')
print -dpng -r1200 tred70.png
end