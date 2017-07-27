close all
addpath(genpath('../../ipic3d_toolbox'));
dir='/data1/gianni/paraview/tred60/';

for cycle=18000:1000:18000
ncyc=num2str(cycle)
leggo=0;
if(leggo==1)

[V,Bx,By,Bz,dx,dy,dz]=read_vtk_3d([dir 'B_xyz_cycle' ncyc '.vtk'],0);
[V,Ex,Ey,Ez,dx,dy,dz]=read_vtk_3d([dir 'E_xyz_cycle' ncyc '.vtk'],0);
[Nx, Ny, Nz] =size(Bx)

%[Bx,By,Bz,Nx,Ny,Nz]=read_binVTK_vector(dir,'B',cycle);
%[Ex,Ey,Ez,Nx,Ny,Nz]=read_binVTK_vector(dir,'E',cycle);
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

%spectrum = fft(Ez.*By-Bz.*Ey,[],3);
spectrum = fft(Bz,[],3);

modes=(squeeze(sum(sum(spectrum.*conj(spectrum),1),2)));
figure(1000)
bar(0:39,(modes(1:40)))
set(gca,'Yscale','log')

for i=1:20
h=figure(i)
subplot(2,1,1)
imagesc(abs(spectrum(:,:,i))')
ylabel('y/d_i')
title(['E_{||}(x,y) m_z=' num2str(i-1) ])
colormap default
colorbar
subplot(2,1,2)
imagesc(angle(spectrum(:,:,i))'.*abs(spectrum(:,:,i))')
load cm_new
colormap(cm_kbwrk)
xlabel('x/d_i')
ylabel('y/d_i')
title(['E_{||}(x,y) m_z=' num2str(i-1) ])
colorbar
%caxis([-2 2]*1e-3)
end
%saveas(h,'tred70.fig')
%print -dpng -r1200 tred70.png
end