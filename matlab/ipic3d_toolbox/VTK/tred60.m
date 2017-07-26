close all
addpath('/home/gianni/matlab2/matlab-parsek');
global contours

contours = 1;

B0=.0097
leggi = 0;
for cyc= 17000:1000:17000
ncyc=num2str(cyc)
if(leggi)
dir='/shared02/gianni/tred60/data/';
[V,Bx,By,Bz,dx,dy,dz]=read_vtk_3d([dir 'B_xyz_cycle' ncyc '.vtk'],0);

[nx,ny,nz]=size(Bx);
[xx yy zz]=ndgrid(1:nx,1:nx, 1:nz);

[V,Ex,Ey,Ez,dx,dy,dz]=read_vtk_3d([dir 'E_xyz_cycle' ncyc '.vtk'],0);

%numvar=2;
%[Rhoe,Rhoi]=read_vtk_multiscalar_3d([dir 'rho_xyz_cycle' ncyc '.vtk'],numvar);


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
ez=squeeze(Ez(:,iy,:))/B0;
rho=squeeze(Rhoi(:,iy,:))*4*pi;

figure('position',[10 10 1000 500])
h(1)=subplot_tight(3,1,1,[0.02,0.01])
pcolor(x,z,by)
shading interp
axis equal
axis tight
g=colorbar
ylabel(g,'B_y/B_0','fontsize',[14])
%ylabel(g,'\rho_i','fontsize',[14])
set(gca,'xtick',[],'xticklabel',[])
set(gca,'fontsize',[14])
ylabel('z/d_i','fontsize',[14])
xlim([0 40])
caxis([-0.6 0.6])
h(2)=subplot_tight(3,1,2,[0.02, 0.01])
pcolor(x,z,ez)
shading interp
axis equal
axis tight
g=colorbar
ylabel(g,'E_z/B_0','fontsize',[14])
set(gca,'fontsize',[14])
xlabel('x/d_i','fontsize',[14])
ylabel('z/d_i','fontsize',[14])
caxis([-0.02 0.02])
xlim([0 40])
linkaxes(h')
set(gcf, 'Renderer', 'Painters');
print('-dpng',['instab_60_' ncyc '.png'])
%close all
end

%ay=vecpot(xx,yy,Bx,By);


%disp('f1')
%h=figure(1);
%set(h,'Position' , [5 5 560 420]);
%coplot(xx,yy,EJe+dEJe+EJi+dEJi,ay,'x/d_i','y/d_i','EJTOT',lmt,[-2e-8 2e-8])
%caxis([-2 2]*1e-8)
%close all



