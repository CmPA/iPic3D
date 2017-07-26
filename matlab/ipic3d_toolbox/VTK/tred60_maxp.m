addpath('/home/gianni/matlab2/matlab-parsek');
global contours
close all

contours = 1;

leggi = 1;
%for cyc= 0:1000:20000
for cyc= 15000:15000
ncyc=num2str(cyc)

dir='/shared02/gianni/tred60/data/'
filename=[dir 'settings.hdf'];

B0 = double(hdf5read(filename,'/collective/Bx0'));
vthe = double(hdf5read(filename,'/collective/species_0/uth'));
vthi = double(hdf5read(filename,'/collective/species_1/uth'));
n0=1/4/pi;
mratio=256;
vthe=.045;


wci=B0;
wpi=1*sqrt(n0);
wce=wci*mratio;
wlh=1/sqrt(1/wce/wci+1/wpi^2);

rhoi=vthi/wci
rhoe=vthe/wce


if(leggi)
[V,Bx,By,Bz,dx,dy,dz]=read_vtk_3d([dir 'B_xyz_cycle' ncyc '.vtk'],0);

[nx,ny,nz]=size(Bx);
[xx yy zz]=ndgrid(1:nx,1:nx, 1:nz);

[V,Ex,Ey,Ez,dx,dy,dz]=read_vtk_3d([dir 'E_xyz_cycle' ncyc '.vtk'],0);

numvar=2;
[Rhoe,Rhoi]=read_vtk_multiscalar_3d([dir 'rho_xyz_cycle' ncyc '.vtk'],numvar);


[V,Jx,Jy,Jz,dx,dy,dz]=read_vtk_3d([dir 'Je_xyz_cycle' ncyc '.vtk'],0);
%numvar=10;
%[Pxx,Pxy,Pxz,Pyy,Pyz,Pzz,Ppar,Pper1,Pper2,Peps]=read_vtk_multiscalar_3d([dir 'Pe_xyz_cycle' ncyc '.vtk'],numvar);

% N.B. x and y are exhanged, damned matlab

lmt=round([nx/8 nx-nx/8 1 ny]);

xx=xx*dx;
yy=yy*dy;
zz=zz*dz;
end

y=1:ny;
for i=1:nx
for k=1:nz
%w=Ppar(i,:,k);
w=Jz(i,:,k).^2;
ymax(i,k)=round(sum(w.*y)./sum(w));
%ymax(i,k)=round(sum(y./V(i,:,k).^2)./sum(1./V(i,:,k).^2));
[dum j] = min((abs(Bx(i,:,k).^2)));
ymax(i,k)=j;

Bymax(i,k)=By(i,round(ymax(i,k)),k)/B0;
Rmax(i,k)=Rhoi(i,round(ymax(i,k)),k)*4*pi;

Ezmax(i,k)=Ez(i,round(ymax(i,k)),k)/B0;

end
end



iy=round(ny/2);

x=squeeze(xx(:,iy,:));
z=squeeze(zz(:,iy,:));
by=squeeze(By(:,iy,:))/B0;
rhoi=squeeze(Rhoi(:,iy,:))/B0;
ez=squeeze(Ez(:,iy,:))/B0;

figure('position',[10 10 500 500])
%h(1)=subplot_tight(4,1,1,[0.02,0.01])
h(1)=subplot(3,1,1)
pcolor(x,z,Rmax)
shading interp
%axis equal
axis tight
g=colorbar
ylabel(g,'n_i/n_0','fontsize',[14])
set(gca,'xtick',[],'xticklabel',[])
set(gca,'fontsize',[14])
ylabel('z/d_i','fontsize',[14])
xlim([20 40])
caxis([0.2 2.2])

%h(2)=subplot_tight(4,1,2,[0.02, 0.01])
h(2)=subplot(3,1,2)
pcolor(x,z,Ezmax)
shading interp
%axis equal
axis tight
g=colorbar
ylabel(g,'E_z/B_0','fontsize',[14])
set(gca,'fontsize',[14])
set(gca,'xtick',[],'xticklabel',[])
ylabel('z/d_i','fontsize',[14])
xlim([20 40])
caxis([-0.02 0.02])

%h(3)=subplot_tight(4,1,3,[0.02, 0.01])
h(3)=subplot(3,1,3)

for i=1:nx
[f, spec] = analisi_fourier(Ezmax(i,:),dz);
nk=max(size(spec));
if(i==1)
SS=zeros(nx,nk);
end
SS(i,1:nk)=spec;
end

pcolor(x(:,1),f*2*pi,log10(abs(SS')))
                         shading interp
                         axis tight
                         %axis square
                         
                         set(gca,'fontsize',[14])
                         %set(gca,'YAxislocation','right')
                         xlabel('x/d_i','fontsize',[14])
                         ylabel('k d_i','fontsize',[14])
                         g=colorbar
                         ylabel(g,'log_{10}(E_{z,FFT})','fontsize',[14])
                         xlim([20 40])
                               ylim([0 20])
                               caxis([-7 -2])

linkaxes(h','x')
%set(gcf, 'Renderer', 'zbuffer');
print('-dpng',['instab_60maxp_' ncyc '.png'])
% close all
end
