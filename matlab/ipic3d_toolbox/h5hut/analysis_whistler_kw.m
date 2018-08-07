close all
clear all
addpath(genpath('../../ipic3d_toolbox'));

folder_name = '/Users/giovannilapenta/Dropbox/Science/codes/build_office/data_full'
namefile = 'Maxwellian';



filename=[folder_name '/settings.hdf'];
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
Bf=zeros(Nx,Nt/Nstep+1);
it=0;

dx=Lx/Nx;
dy=Ly/Ny;
dz=Lz/Nz;


[X Y] = meshgrid(0:dx:Lx-dx,0:dy:Ly-dy);
sigma=Ly/5
Kernel1=exp(-((X-Lx/2).^2+(Y-Ly/2).^2)/sigma.^2)';

Bf=[];Ef=[];
for cycle= 0:Nstep:Nt
   

    itn=sprintf('%06.0f',cycle);
        
    fn=[folder_name,'/',[namefile '-Fields'],'_',itn,'.h5']

    hinfo=hdf5info(fn);
    Nx= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(1);
    Ny= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(2);
    Nz= hinfo.GroupHierarchy.Groups.Groups.Groups(3).Datasets(1).Dims(3);
    % uncomment this for a list of varibales available
    %hinfo.GroupHierarchy.Groups.Groups.Groups(:).Name
    
    
old=0;
    if(old)
    bx = hdf5read(fn,'/Step#0/Block/Btx/0/');
    by = hdf5read(fn,'/Step#0/Block/Bty/0/');
    bz = hdf5read(fn,'/Step#0/Block/Btz/0/');
    else
    bx = hdf5read(fn,'/Step#0/Block/Bx/0/');
    by = hdf5read(fn,'/Step#0/Block/By/0/');
    bz = hdf5read(fn,'/Step#0/Block/Bz/0/');
    bx_ext = hdf5read(fn,'/Step#0/Block/Bx_ext/0/');
    by_ext = hdf5read(fn,'/Step#0/Block/By_ext/0/');
    bz_ext = hdf5read(fn,'/Step#0/Block/Bz_ext/0/');
    bx=bx+bx_ext;
    by=by+by_ext;
    bz=bz+bz_ext;
    end
    
  
    ex = hdf5read(fn,'/Step#0/Block/Ex/0/');
    ey = hdf5read(fn,'/Step#0/Block/Ey/0/');
    ez = hdf5read(fn,'/Step#0/Block/Ez/0/');

  
    
    b = sqrt (bx.^2 +by.^2 + bz.^2);
    
   
     epar=dot(ex,ey,ez,bx,by,bz)./b;
Bf=[Bf squeeze(by(:,1,1))];%.*Kernel1;
Ef=[Ef squeeze(ex(:,1,1))];%.*Kernel1;
end


VA=B0;
wce=B0*abs(qom);
wpe=sqrt(abs(qom));
di=1.0;
de=di/sqrt(abs(qom));


kx=pi/dx*linspace(-1,1,Nx);
ky=pi/dy*linspace(-1,1,Ny);
omega=pi/(Dt*Nstep)*linspace(-1,1,Nt/Nstep+1);

%Bfft=fft(fft(fft(Bf,[],1),[],2),[],3);
Bfft=fftshift(fft2(Bf));
Efft=fftshift(fft2(Ef));

figure(1)
vec=abs(Bfft);
%vec=vec+flipud(vec);
%vec=vec+fliplr(vec);
imagesc(kx*de,omega/wpe,log10(vec)')
% xlim([-1 1])
% ylim([-1 1])
hold on 
%plot(ky*de,(pi*di*VA*(ky).^2./(1+(ky*de).^2))/wce,'w','linewidth',2);
%plot(ky*de,((ky*de).^2./(1+(ky*de).^2)),'w--','linewidth',2);
%caxis([-6 2])
xlabel('k_{||}d_e')
ylabel('\omega/\Omega_{pe}')
colormap hsv
caxis(max(max(log10(vec)))+[-3 0])
colorbar
caxis([-3 2.5])
colormap jet
print('-dpng','-r300',[folder_name 'Bkpar_kperp0_w'  '.png'])


figure(2)
vec=abs(Efft);
%vec=vec+flipud(vec);
%vec=vec+fliplr(vec);
imagesc(kx*de,omega/wpe,log10(vec)')
% xlim([-1 1])
% ylim([-1 1])
hold on 
%plot(ky*de,(pi*di*VA*(ky).^2./(1+(ky*de).^2))/wce,'w','linewidth',2);
%plot(ky*de,((ky*de).^2./(1+(ky*de).^2)),'w--','linewidth',2);
%caxis([-6 2])
xlabel('k_{||}d_e')
ylabel('\omega/\Omega_{pe}')
colormap hsv
caxis(max(max(log10(vec)))+[-3 0])
colorbar
caxis([-3 2.5])
colormap jet
print('-dpng','-r300',[folder_name 'Ekpar_kperp0_w'  '.png'])
return
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