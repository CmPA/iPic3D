maxp_common


iz = Nz-round(Nz*(max(Zgsmrange)-Ygsm)/(max(Zgsmrange)-min(Zgsmrange)));
Ygsm=gsmz2y(Lz-zc(iz));



for cycle=Ncyc_ini:1000:Ncyc_max

ncycle=num2str(cycle)
ncycle1=num2str(cycle,'%06d')


global blowup contours
blowup=0;
contours=1;


file=[dir 'B_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
Bx=reshape(V,Nx,Ny,Nz);


file=[dir 'B_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
By=reshape(V,Nx,Ny,Nz);


file=[dir 'B_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
Bz=reshape(V,Nx,Ny,Nz);


file=[dir 'rho_1_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
rhoi=reshape(V,Nx,Ny,Nz);


for ix=1:Nx
close all

%ix=round(Nx/3);
nix=num2str(ix,'%04d');


BxYZ=squeeze(Bx(ix,:,:));

ByYZ=squeeze(By(ix,:,:));

BzYZ=squeeze(Bz(ix,:,:));

rhoiYZ=squeeze(rhoi(ix,:,:));

BYZ=sqrt(BxYZ.^2+ByYZ.^2+BzYZ.^2);

n0=4*pi*(max(rhoiYZ(:))+min(rhoiYZ(:)))/2;
n0=min(rhoiYZ(:))*4*pi;
di=1/sqrt(n0)
de=di/sqrt(256)
ri=vthi/mean(BYZ(:))
re=-vthe/mean(BYZ(:))/qom
rlh=sqrt(ri*re)
B0=5% This is the harmonic mean of B
re=-vthe./(B0/code_B)/qom
ri=vthi/mean(B0/code_B)

Fs = 1/dy;                    % Sampling frequency
T = 1/Fs;                     % Sample length
ky = Fs*pi*linspace(0,1,Ny/2+1);


Fs = 1/dz;                    % Sampling frequency
T = 1/Fs;                     % Sample length
kz = Fs*pi*linspace(0,1,Nz/2+1);

k=ky*sqrt(2);
fBx = fft2(BxYZ.^2);
fBx=fBx(1:Ny/2+1,1:Nz/2+1);
%pcolor(ky,kz,log(abs(fBx)))
%shading interp
spec=diag(abs(fBx));
loglog(k,spec,'k','LineWidth',2)
hold on
loglog(k,spec(2).*k.^(-5/3)./k(2).^(-5/3),'y',k,spec(end).*k.^(-5)./k(end).^(-5),'r')
span=[min(spec) max(spec)];
loglog(1./[di di],span,'g--',1./[de de],span,'g',1./[rlh rlh],span,'m',1./[re re],span,'b',1./[ri ri],span,'b--')
xlabel('kd_{i0}','fontsize',[14])
ylabel('FFT(B_x^2)','fontsize',[14])
legend('spectrum(Bx)','\gamma=-5/3','\gamma=-5','d_i^{-1}','d_e^{-1}','\rho_{LH,loc}^{-1}','\rho_{e0}^{-1}','\rho_{i0}^{-1}','location','southwest')
set(gca,'fontsize',[14])
ylim([1e-7 1e-1])
title(['X(gsm)/R_E= ' num2str(gsmx(xc(ix)))],'fontsize',[14])
set(gcf,'Renderer','zbuffer');
print('-dpng',['k_spectrum_' nix '.png'])
end

end
