close all
clear all
addpath(genpath('../../ipic3d_toolbox'));
dir='/data1/gianni/iPic3D-github/build/data/';
%dir='/data2/gianni/whistler/whistler25/';

filename=[dir 'settings.hdf'];
Dt=double(hdf5read(filename,'/collective/Dt'));
B0=double(hdf5read(filename,'collective/By0'));
Lx=double(hdf5read(filename,'/collective/Lx'));
Ly=double(hdf5read(filename,'/collective/Ly'));
Lz=double(hdf5read(filename,'/collective/Lz'));
Nx=double(hdf5read(filename,'/collective/Nxc')); 
Ny=double(hdf5read(filename,'/collective/Nyc'));
Nz=double(hdf5read(filename,'/collective/Nzc'));
XLEN=double(hdf5read(filename,'/topology/XLEN'));
YLEN=double(hdf5read(filename,'/topology/YLEN'));
ZLEN=double(hdf5read(filename,'/topology/ZLEN'));
qom=double(hdf5read(filename,'/collective/species_0/qom'));
vthi=double(hdf5read(filename,'/collective/species_1/uth'));
vthe=double(hdf5read(filename,'/collective/species_0/uth'));
Nprocs=hdf5read(filename,'/topology/Nprocs');

%for cycle=0:20:3760

Nt=400;
Nstep=2;
%Nt=398
%Nstep=2
Bf=zeros(Nx,Ny,Nt/Nstep+1);
it=0;

dx=Lx/Nx;
dy=Ly/Ny;
dz=Lz/Nz;


[X Y] = meshgrid(0:dx:Lx-dx,0:dy:Ly-dy);
sigma=Ly/5
Kernel1=exp(-((X-Lx/2).^2+(Y-Ly/2).^2)/sigma.^2)';

for cycle= 0:Nstep:Nt
    it=it+1;

[Bx,By,Bz,Nx,Ny,Nz]=read_binVTK_vector(dir,'B',cycle);
[Ex,Ey,Ez,Nx,Ny,Nz]=read_binVTK_vector(dir,'E',cycle);
Bf(:,:,it)=Bz;%.*Kernel1;
%[Ex,Ey,Ez,Nx,Ny,Nz]=read_binVTK_vector(dir,'E',cycle);
%[Jex,Jey,Jez,Nx,Ny,Nz]=read_binVTK_vector(dir,'Je',cycle);
%[Jix,Jiy,Jiz,Nx,Ny,Nz]=read_binVTK_vector(dir,'Ji',cycle);

%[Az,Nx,Ny,Nz]=read_binVTK_scalar(dir,'Az',cycle);
%[rhoe,rhoi,Nx,Ny,Nz]=read_binVTK_multiscalar(dir,'rho',cycle);
%[Pixx,Pixy,Pixz,Piyy,Piyz,Pizz,Pipar,Piper1,Piper2,Pieps,Nx,Ny,Nz] = read_binVTK_pressure(dir,'Pi',cycle);
%[Pexx,Pexy,Pexz,Peyy,Peyz,Pezz,Pepar,Peper1,Peper2,Peeps,Nx,Ny,Nz] = read_binVTK_pressure(dir,'Pe',cycle);

%B=sqrt(Bx.*Bx+By.*By+Bz.*Bz);
%Epar=(Ex.*Bx+Ey.*By+Ez.*Bz)./B;

%Te=(Pexx+Peyy+Pezz)./(-rhoe);
%Ti=(Pixx+Piyy+Pizz)./rhoi;
end

VA=B0;
wce=B0*abs(qom);
di=1.0;
de=di/sqrt(abs(qom));


kx=pi/dx*linspace(-1,1,Nx);
ky=pi/dy*linspace(-1,1,Ny);
omega=pi/(Dt*Nstep)*linspace(-1,1,Nt/Nstep+1);

%Bfft=fft(fft(fft(Bf,[],1),[],2),[],3);
Bfft=fftshift(fftn(Bf));

figure(1)
vec=squeeze(sum(abs(Bfft(end/2-1:end/2+1,:,:)).^2,1))';
vec=vec+flipud(vec);
vec=vec+fliplr(vec);
imagesc(kx*de,omega/wce,log10(vec))
% xlim([-1 1])
% ylim([-1 1])
hold on 
%plot(ky*de,(pi*di*VA*(ky).^2./(1+(ky*de).^2))/wce,'w','linewidth',2);
plot(ky*de,((ky*de).^2./(1+(ky*de).^2)),'w--','linewidth',2);
%caxis([-6 2])
xlabel('k_{||}d_e')
ylabel('\omega/\Omega_{ce}')
colormap hsv
caxis(max(max(log10(vec)))+[-3 0])
colorbar
caxis([-3 2.5])
colormap jet
print('-dpng','-r300',[dir 'kpar_kperp0_w'  '.png'])

% figure(100)
% vec=squeeze(sum(abs(Bfft(end/4-10:end/4+10,:,:)).^2,1))';
% vec=vec+flipud(vec);
% vec=vec+fliplr(vec);
% imagesc(kx*de,omega/wce,log10(vec))
% % xlim([-1 1])
% % ylim([-1 1])
% %caxis([-6 2])
% hold on 
% %plot(ky*de,(pi*di*VA*(ky).^2./(1+(ky*de).^2))/wce,'w','linewidth',2);
% kxx=kx(end/4)
% plot(ky*de,((ky*de).*sqrt((kxx*de).^2+(kxx*de).^2)./(1+(ky*de).^2+(kxx*de).^2)),'w--','linewidth',2);
% caxis(max(max(log10(vec)))+[-3 0])
% xlabel('k_{||}d_e')
% ylabel('\omega/\Omega_{ce}')
% colormap hsv
% colorbar
% print('-dpng','-r300',[dir 'kpar2_kperp_w'  '.png'])

figure(2)
w=log10(squeeze(sum(abs(Bfft).^2,1))')
imagesc(ky*de,omega/wce,w)

%xlim([-1 1])
%ylim([-1 1])
hold on 
%plot(ky*de,(pi*di*VA*(ky).^2./(1+(ky*de).^2))/wce,'w','linewidth',2);
plot(ky*de,((ky*de).^2./(1+(ky*de).^2)),'w--','linewidth',2);
caxis(max(max(w))+[-3 0])
%caxis(4.0+[-3 0])

xlabel('k_{||}d_e')
ylabel('\omega/\Omega_{ce}')
caxis([-3 2.5])
colormap jet
colorbar
print('-dpng','-r300',[dir 'kpar_int_kperp_w'  '.png'])


figure(3)
imagesc(kx*de,omega/wce,log10(squeeze(sum(abs(Bfft).^2,2))'))
% xlim([-1 1])
% ylim([-1 1])
xlabel('k_{\perp}d_e')
ylabel('\omega/\Omega_{ce}')
caxis(max(max(log10(squeeze(sum(abs(Bfft).^2,2))')))+[-3 0])
caxis([-3 2.5])
colormap jet
colorbar
print('-dpng','-r300',[dir 'kper_w'  '.png'])

%close all
figure(4)

w=log10(squeeze(sum(abs(Bfft(end/2-1:end/2+1,:,:)).^2,1))');

wwhis=(ky*de).^2./(1+(ky*de).^2);
[kk,ww]=meshgrid(ky.*de,omega/wce);
pcolor(kk,ww./kk,w); shading interp; caxis(max(max(w))+[-3 0]);
caxis([-3 2.5])
colormap jet
colorbar
xlim([-7 7])
ylim([-2 2])
hold on
plot(ky*de,wwhis./(ky*de),'w',ky*de,-wwhis./(ky*de),'w',ky*de,(1-wwhis)./(ky*de),'k',ky*de,(1+wwhis)./(ky*de),'k','LineWidth',2)
xlabel('k_{||}d_e')
ylabel('\omega/(k_{||} \omega_{ce} d_e)')
figure(4)
print('-dpng','-r300',[dir 'komega_nulkperp.png'])



figure(5)
w=log10(squeeze(sum(abs(Bfft).^2,1))')


wwhis=(ky*de).^2./(1+(ky*de).^2);
[kk,ww]=meshgrid(ky.*de,omega/wce);
pcolor(kk,ww./kk,w); shading interp; caxis(max(max(w))+[-3 0]);
caxis([-3 2.5])
colormap jet
colorbar
xlim([-7 7])
ylim([-2 2])
hold on
plot(ky*de,wwhis./(ky*de),'w',ky*de,-wwhis./(ky*de),'w',ky*de,(1-wwhis)./(ky*de),'k',ky*de,(1+wwhis)./(ky*de),'k','LineWidth',2)
figure(5)
xlabel('k_{||}d_e')
ylabel('\omega/(k_{||} \omega_{ce} d_e)')
print('-dpng','-r300',[dir 'komega_intkperp.png'])