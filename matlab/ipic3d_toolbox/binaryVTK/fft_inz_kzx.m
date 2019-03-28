close all
addpath(genpath('../../ipic3d_toolbox'));


sim_name='tred81'
switch sim_name
case 'tred77'
TRED77;
case_name='GEM';
cycle = 15000;
zcode = Lz/2;
case 'tred81'
tred81;
case_name='GEM';
cycle = 18000;
zcode = Lz/2;
    case 'tred82'
tred82;
case_name='GEM';
cycle = 18000;
zcode = Lz/2;    
    case 'tred85'
tred85;
case_name='GEM';
cycle = 18000
zcode = Lz/2;
case 'AH'
generic;
case_name='AH';
cycle =4000;
zcode = Lz/2;
case 'HRmaha3D3'
HRmaha3D3;
    case_name='GEM';
dir='/data1/gianni/HRmaha3D3/h5/'; cycle= 80002; ncycle = num2str(cycle,'%06d');
cycle = 80002;  % for h5
%cycle = 80000  % for vtk binary
% for HRmaha3D1:
time=60*(cycle/75000.0*Dt/.125); %*4 %times four to correct for change in dt between 2D and 3D
% for HRmaha3D1.v2:
% time=60*(cycle/75000.0) *2 %times two to correct for change in dt between 2D and 3D
%ADD initial time of the RUN
time=time+initial_time; %(03*60+48)*60
case '7feb09'
FEB09;
cycle=18000
case_name='MHDUCLA'
%cycle = 80000  % for vtk binary
% for HRmaha3D1:
time=60*(cycle/75000.0*Dt/.125); %*4 %times four to correct for change in dt between 2D and 3D
% for HRmaha3D1.v2:
% time=60*(cycle/75000.0) *2 %times two to correct for change in dt between 2D and 3D
%ADD initial time of the RUN
time=time+initial_time; %(03*60+48)*60
otherwise
print('no recognised case selected')
end

% Prepare string
ntime = num2str(cycle,'%06d');
ncycle = num2str(cycle,'%06d');

must_read=true;
leggo='h5';
if(must_read)
import_h5_binvtk
end

% for cycle=18000:1000:18000
% ncyc=num2str(cycle)
% leggo=0;
% if(leggo==1)
% 
% [V,Bx,By,Bz,dx,dy,dz]=read_vtk_3d([dir 'B_xyz_cycle' ncyc '.vtk'],0);
% [V,Ex,Ey,Ez,dx,dy,dz]=read_vtk_3d([dir 'E_xyz_cycle' ncyc '.vtk'],0);
% [Nx, Ny, Nz] =size(Bx)
% 
% %[Bx,By,Bz,Nx,Ny,Nz]=read_binVTK_vector(dir,'B',cycle);
% %[Ex,Ey,Ez,Nx,Ny,Nz]=read_binVTK_vector(dir,'E',cycle);
% % [Jex,Jey,Jez,Nx,Ny,Nz]=read_binVTK_vector(dir,'Je',cycle);
% % [Jix,Jiy,Jiz,Nx,Ny,Nz]=read_binVTK_vector(dir,'Ji',cycle);
% % 
% % [Az,Nx,Ny,Nz]=read_binVTK_scalar(dir,'Az',cycle);
% % [rhoe,rhoi,Nx,Ny,Nz]=read_binVTK_multiscalar(dir,'rho',cycle);
% % [Pixx,Pixy,Pixz,Piyy,Piyz,Pizz,Pipar,Piper1,Piper2,Pieps,Nx,Ny,Nz] = read_binVTK_pressure(dir,'Pi',cycle);
% % [Pexx,Pexy,Pexz,Peyy,Peyz,Pezz,Pepar,Peper1,Peper2,Peeps,Nx,Ny,Nz] = read_binVTK_pressure(dir,'Pe',cycle);
% % 
%  B=sqrt(Bx.*Bx+By.*By+Bz.*Bz);
%  Epar=(Ex.*Bx+Ey.*By+Ez.*Bz)./B;
% % 
% % Te=(Pexx+Peyy+Pezz)./(-rhoe);
% % Ti=(Pixx+Piyy+Pizz)./rhoi;
% end

% Lx=40;Ly=20;Lz=10;
% 
% dx=Lx/Nx;
% dy=Ly/Ny;
% dz=Lz/Nz;

[X Y] = meshgrid(dx/2:dx:Lx-dx/2,dy/2:dy:Ly-dy/2);

[X Y] = meshgrid(0:dx:Lx-dx,0:dy:Ly-dy);

[Sx, Sy, Sz] = cross_prod(Ex, Ey, Ez, Bx/4/pi, By/4/pi, Bz/4/pi);
S=sqrt(Sx.^2+Sy.^2+Sz.^2);
[x,y,z]=meshgrid(0:dx:Lx,0:dy:Ly,0:dz:Lz);
radius=0.01;
divS = compute_div(x,y,z,Sx,Sy,Sz,radius, 0.);

xc=Lx-linspace(0, Lx, Nx);
yc=linspace(0, Ly, Ny);
AAz=zeros(size(Bx));
for kr=1:Nz
AAz(:,:,kr)=vecpot(xc,yc,Bx(:,:,kr),By(:,:,kr));
AAz(:,:,kr)=AAz(:,:,kr)-AAz(round(Nx/2),round(Ny/2),kr);
end
figure
imagesc(mean(AAz,3)')
colorbar

[X,Y,Z]=ndgrid(1:Nx,1:Ny,1:Nz);
figure
plot3(AAz(:),Y(:),S(:),'.')
figure
plot(AAz(:),S(:),'.')
figure
ndiv=100;
[totnum,nbinu,xrange,urange]=spaziofasi2(AAz(:),S(:),ones(Nx*Ny*Nz,1),0,min(AAz(:)),max(AAz(:)),min(S(:)),max(S(:))/10,ndiv);
imagesc(xrange,urange,log10(nbinu))
close all

skippa=false
if(skippa)
work_horse_y(Sx,'Sx',Lx,Ly);
work_horse_y(Sy,'Sy',Lx,Ly);
work_horse_y(Sz,'Sz',Lx,Ly);
work_horse_y(S,'S',Lx,Ly);
work_horse_y(divS,'divS',Lx,Ly);
work_horse_y(Ex,'Ex',Lx,Ly);
work_horse_y(Ey,'Ey',Lx,Ly);
work_horse_y(Ez,'Ez',Lx,Ly);
work_horse_y(Bx,'Bx',Lx,Ly);
work_horse_y(By,'By',Lx,Ly);
work_horse_y(Bz,'Bz',Lx,Ly);
end
work_horse_y(divS,'divS',Lx,Ly,Lz);
        
        Epx = Ex + (Jey.*Bz - Jez.*By)./rhoe;
        Epy = Ey + (Jez.*Bx - Jex.*Bz)./rhoe;
        Epz = Ez + (Jex.*By - Jey.*Bx)./rhoe;
        
        JdotEp=(Jex+Jix).*Epx + (Jey+Jiy).*Epy + (Jez+Jiz).*Epz;
work_horse_y(imgaussfilt3(JdotEp,2.0),'JEp',Lx,Ly,Lz);
        
%saveas(h,'tred70.fig')
%print -dpng -r1200 tred70.png
%end

function [out]=work_horse(S,varname,Lx,Ly)

%spectrum = fft(Ez.*By-Bz.*Ey,[],3);
spectrum = fft(S,[],3);

modes=(squeeze(sum(sum(spectrum.*conj(spectrum),1),2)));

for i=1:20
    nc(i)=round(mean(mean(log10(1e-10+abs(spectrum(:,:,i))'))));
end
nc(2:end)=max(nc(2:end));
figure(1000)
bar(0:39,(modes(1:40)))
set(gca,'Yscale','log')
close all
for i=1:20
h=figure(100+i)
%subplot(2,1,1)
imagesc([0 Lx],[0 Ly],log10(abs(spectrum(:,:,i))'))
ylabel('y/d_i')
title([varname '(x,y) m_z=' num2str(i-1) ])

nc(i)
caxis([nc(i)-2 nc(i)+2])

colormap hsv
colorbar
axis equal 
axis tight
% subplot(2,1,2)
% imagesc(angle(spectrum(:,:,i))')%.*abs(spectrum(:,:,i))')
% load cm_new
% %colormap(cm_kbwrk)
% xlabel('x/d_i')
% ylabel('y/d_i')
% title([varname '(x,y) m_z=' num2str(i-1) ])
% colorbar
print(h,'-dpng', '-r300',['FFTZ_' varname '_m' num2str(i-1,'%03.f') '.png'])
%export_fig(['eFFTZ_' varname '_m' num2str(i-1,'%03.f') '.png'],'-png')
close(100+i)
%caxis([-2 2]*1e-3)
end
out=1;
end

function [out]=work_horse_y(S,varname,Lx,Ly,Lz)

%spectrum = fft(Ez.*By-Bz.*Ey,[],3);
spectrum = fft(S,[],3);
a=(abs(squeeze(spectrum(:,96,:)))');
a=a(2:40,:);
V=a*(diag(1./sum(a)))
imagesc([0 Lx], [2-1 40-1]*2*pi/Lz, log(a*(diag(1./sum(a)))))
load gist_ncar
colormap(gist_ncar)
set(gca,'fontsize',[14])
xlabel('x/d_i','fontsize',[14])
ylabel('k_zd_i','fontsize',[14])


title(varname,'fontsize',[14])

axis image
axis xy
colorbar
caxis([-10 0])
% subplot(2,1,2)
% imagesc(angle(spectrum(:,:,i))')%.*abs(spectrum(:,:,i))')
% load cm_new
% %colormap(cm_kbwrk)
% xlabel('x/d_i')
% ylabel('y/d_i')
% title([varname '(x,y) m_z=' num2str(i-1) ])
% colorbar
print('-dpng', '-r300',['FFTZ_Y_' varname '_m' num2str(i-1,'%03.f') '.png'])

semilogy(linspace(1,39,39)*2*pi/Lz, mean(V,2))
ylabel(['FFT - ' varname] ,'fontsize',[14])
xlabel('k_zd_i','fontsize',[14])
print('-dpng', '-r300',['FFTZ_YX_' varname '_m' num2str(i-1,'%03.f') '.png'])
ciplot(min(log10(V')),max(log10(V')),linspace(1,39,39)*2*pi/Lz,'yellow')
ylabel(['FFT - ' varname] ,'fontsize',[14])
xlabel('k_zd_i','fontsize',[14])
print('-dpng', '-r300',['FFTZ_YXband_' varname '_m' num2str(i-1,'%03.f') '.png'])
%export_fig(['eFFTZ_' varname '_m' num2str(i-1,'%03.f') '.png'],'-png')
%caxis([-2 2]*1e-3)
end
