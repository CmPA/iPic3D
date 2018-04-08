function [] = spettro_alongz(x,y,z,Ex,Ey,Ez, Nx, Ny, Nz,Lx,Ly,Lz,vthi,B0)


close all
h1=figure(1)
set(h1,'Position',[40 378 1040 420])

cmap=load('rainbow_cm.mat');
cmap=cmap.cmap;


[cEx,cEy,cEz]=compute_curl(x,y,z,Ex,Ey,Ez);

fEx = fft(Ex(:,:,1:end-1),[],3);fEx=squeeze(mean(abs(fEx(:,:,1:round(Nz/2)+1)),2));
fEy = fft(Ey(:,:,1:end-1),[],3);fEy=squeeze(mean(abs(fEy(:,:,1:round(Nz/2)+1)),2));
fEz = fft(Ez(:,:,1:end-1),[],3);fEz=squeeze(mean(abs(fEz(:,:,1:round(Nz/2)+1)),2));

max_val=-7;%log(max(max(fEx(10:end-10,:))));
min_val=-10;

x1d=10:Nx-10; x1d=x1d/Nx*Lx;

mz=0:round(Nz/2); kz=mz*2*pi/Lz;

rhoi=vthi/B0
rhoe=rhoi/sqrt(256)
rholh=sqrt(rhoi*rhoe)

% subplot(2,3,1)
% pcolor(x1d,kz,log(fEx(10:end-10,:)'));colorbar;shading interp 
% caxis([-10 max_val])
% colormap(jet)
% xlabel('x/d_i')
% ylabel('k_z d_i')
% ylim([0 20])

grafo(fEx,1,'Spectrum Ex',cmap)
grafo(fEz,2,'Spectrum Ey',cmap)
grafo(fEy,3,'Spectrum Ez',cmap)



fEx = fft(cEx(:,:,1:end-1),[],3);fEx=squeeze(mean(abs(fEx(:,:,1:round(Nz/2)+1)),2));
fEy = fft(cEy(:,:,1:end-1),[],3);fEy=squeeze(mean(abs(fEy(:,:,1:round(Nz/2)+1)),2));
fEz = fft(cEz(:,:,1:end-1),[],3);fEz=squeeze(mean(abs(fEz(:,:,1:round(Nz/2)+1)),2));

grafo(fEx,4,'Spectrum dBx/dt',cmap)
grafo(fEz,5,'Spectrum dBy/dt',cmap)
grafo(fEy,6,'Spectrum dBz/dt',cmap)

print -dpng figura

function [] = grafo(fEx,pos,varname,cmap)
subplot(2,3,pos)
pcolor(x1d,kz,log(fEx(10:end-10,:)'));colorbar;shading interp 
caxis([-10 max_val])
xlabel('x/d_i')
ylabel('k_y d_i')
ylim([0 20])
hold on 
plot(x1d, ones(size(x1d))/rhoi,'k','linewidth',2)
plot(x1d, ones(size(x1d))/rhoe,'k:','linewidth',2)
plot(x1d, ones(size(x1d))/rholh,'k--','linewidth',2)
title(varname)
colormap(flipud(cmap))
end

end