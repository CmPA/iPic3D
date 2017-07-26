addpath('/home/gianni/matlab2/matlab-parsek');
global contours
close all

contours = 1;

dir0= '/shared/gianni/tred54/'
filename=[dir0 'settings.hdf'];
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


leggi = 0;
for cyc= 11000:1000:11000
ncyc=num2str(cyc)
if(leggi)

dir= [dir0 'vtk-new/']
if (cyc>12000)
dir='/shared/gianni/tred54.2/vtk/'
end

[V,Bx,By,Bz,dx,dy,dz]=read_vtk_3d([dir 'fields_B_cycle' ncyc '.vtk'],1);

[nx,ny,nz]=size(Bx);
[xx yy zz]=ndgrid(1:nx,1:nx, 1:nz);

[V,Ex,Ey,Ez,dx,dy,dz]=read_vtk_3d([dir 'fields_E_cycle' ncyc '.vtk'],1);

[Rhoi,Dx,Dy,Dz,dx,dy,dz]=read_vtk_3d([dir 'rho_tot_ions_cycle' ncyc '.vtk'],1);

% N.B. x and y are exhanged, damned matlab

lmt=round([nx/8 nx-nx/8 1 ny]);

xx=xx*dx;
yy=yy*dy;
zz=zz*dz;
end

iy=round(ny/2);

x=squeeze(xx(:,iy,:));
z=squeeze(zz(:,iy,:));
by=squeeze(By(:,iy,:))/B0;
rhoi=squeeze(Rhoi(:,iy,:));
ez=squeeze(Ez(:,iy,:))/B0;

figure('position',[10 10 250 500])
%h(1)=subplot_tight(4,1,1,[0.02,0.01])
h(1)=subplot(3,1,1)
pcolor(x,z,rhoi)
shading interp
%axis equal
axis tight
%g=colorbar
%ylabel(g,'B_y/B_0','fontsize',[14])
set(gca,'xtick',[],'xticklabel',[])
%set(gca,'YAxislocation','right')
set(gca,'fontsize',[14])
ylabel('z/d_i','fontsize',[14])
xlim([0 10])
caxis([0.2 2.2])

%h(2)=subplot_tight(4,1,2,[0.02, 0.01])
h(2)=subplot(3,1,2)
pcolor(x,z,ez)
shading interp
%axis equal
axis tight
%g=colorbar
%ylabel(g,'E_z/B_0','fontsize',[14])
set(gca,'fontsize',[14])
set(gca,'xtick',[],'xticklabel',[])
%set(gca,'YAxislocation','right')
%xlabel('x/d_i','fontsize',[14])
ylabel('z/d_i','fontsize',[14])
xlim([0 10])
caxis([-0.02 0.02])


%h(3)=subplot_tight(4,1,3,[0.02, 0.01])
h(3)=subplot(3,1,3)

for i=1:nx
[f, spec] = analisi_fourier(ez(i,:),dz);
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
                         %g=colorbar
                         %ylabel(g,'log_{10}(E_{z,FFT})','fontsize',[14])
                         xlim([0 10])
                         ylim([0 20])
                         caxis([-7 -2])
                         
linkaxes(h','x')
set(gcf, 'Renderer', 'zbuffer');

print('-dpng',['instab_' ncyc '.png'])

end
