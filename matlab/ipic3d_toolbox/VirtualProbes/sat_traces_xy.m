addpath 'matlab-parsek'

close all
clear all

vthr=.045;

results_dir='/shared/gianni/tred54.2/';

filename=[results_dir 'settings.hdf']
Lx=hdf5read(filename,'/collective/Lx')
Ly=hdf5read(filename,'/collective/Ly')
Lz=hdf5read(filename,'/collective/Lz')
B0x=hdf5read(filename,'/collective/Bx0')
Dt=hdf5read(filename,'/collective/Dt')
XLEN=hdf5read(filename,'/topology/XLEN')
YLEN=hdf5read(filename,'/topology/YLEN')
ZLEN=hdf5read(filename,'/topology/ZLEN')
mratio=abs(hdf5read(filename,'/collective/species_0/qom'))
Nprocs=hdf5read(filename,'/topology/Nprocs')

% instability in cavities
ipx=XLEN/2-1; ipy=YLEN/2+2*0;
ipx=9; ipy=YLEN/2+1;
%ipx=XLEN/2+3; ipy=YLEN/2;

%ipx=XLEN/2; ipy=YLEN/2;
%ipx=XLEN/2+4; ipy=YLEN/2+1;

%whistler in DF
ipx = XLEN/2 -3
ipy = YLEN/2+1

zzp=[];
EEX=[];
EEY=[];
EEZ=[];
BBX=[];
BBY=[];
BBZ=[];
JPAR=[];
JX=[];
JY=[];
JZ=[];
JXI=[];
JYI=[];
JZI=[];
NE=[];
NI=[];
for ipz=1:ZLEN;
ip=(ipx-1)*YLEN*ZLEN+(ipy-1)*ZLEN+ipz-1;

nome=[results_dir 'VirtualSatelliteTraces' num2str(ip) '.txt']
system(['gunzip ' nome])

fid=fopen(nome);
for i=1:27
x=fscanf(fid,'%f',3); 
xp(i)=x(1);
yp(i)=x(2);
zp(i)=x(3);
end
mean(xp)
mean(yp)
mean(zp)


a=fscanf(fid,'%f',[14 inf])';

fclose(fid);
skip=0;
bx=a(:,1+skip);
by=a(:,2+skip);
bz=a(:,3+skip);
ex=a(:,4+skip);
ey=a(:,5+skip);
ez=a(:,6+skip);
jxe=a(:,7+skip);
jye=a(:,8+skip);
jze=a(:,9+skip);
jxi=a(:,10+skip);
jyi=a(:,11+skip);
jzi=a(:,12+skip);
rhoe=a(:,13+skip)*4*pi;
rhoi=a(:,14+skip)*4*pi;

b=sqrt(bx.*bx+by.*by+bz.*bz);
epar=(ex.*bx+ey.*by+ez.*bz)./b;
jepar=(jxe.*bx+jye.*by+jze.*bz)./b;

[n m]=size(bx);

n0=mean(rhoi(:)-rhoe(:))/2; n0_cavity=.02
n0_cavity=n0
b0=mean(sqrt(bx(:).^2+by(:).^2+bz(:).^2));%b0=B0x;
wci=b0;
wpi=1*sqrt(n0);
wpi_cavity=1*sqrt(n0_cavity);
wce=wci*mratio;
wlh=1/sqrt(1/wce/wci+1/wpi_cavity^2);
%wpi=1 %apparently the plasma oscillations are generated elsewhere where n0=1

n1=floor(n/27)
%n1=14555
%n1=6851
%n1=12530

%n1=14000;

n1=3200

ex=reshape(ex(1:n1*27),27,n1);
ey=reshape(ey(1:n1*27),27,n1);
ez=reshape(ez(1:n1*27),27,n1);
bx=reshape(bx(1:n1*27),27,n1);
by=reshape(by(1:n1*27),27,n1);
bz=reshape(bz(1:n1*27),27,n1);
rhoe=reshape(rhoe(1:n1*27),27,n1);
rhoi=reshape(rhoi(1:n1*27),27,n1);
jepar=reshape(jepar(1:n1*27),27,n1);
jze=reshape(jze(1:n1*27),27,n1);
jye=reshape(jye(1:n1*27),27,n1);
jxe=reshape(jxe(1:n1*27),27,n1);
jzi=reshape(jzi(1:n1*27),27,n1);
jyi=reshape(jyi(1:n1*27),27,n1);
jxi=reshape(jxi(1:n1*27),27,n1);
t=linspace(0,n1,n1);




isatx=1;isaty=1;isatz=1:3;
isat=(isatx-1)*3*3+(isaty-1)*3+isatz;

%[xp(isat);yp(isat);zp(isat)]
zzp=[zzp zp(isat)]
EEX=[EEX ;ex(isat,:)];
EEY=[EEY ;ey(isat,:)];
EEZ=[EEZ ;ez(isat,:)];
BBX=[BBX ;bx(isat,:)];
BBY=[BBY ;by(isat,:)];
BBZ=[BBZ ;bz(isat,:)];

JPAR=[JPAR ;jepar(isat,:)./rhoe(isat,:)];
JZ=[JZ ;jze(isat,:)./rhoe(isat,:)];
JY=[JY ;jye(isat,:)./rhoe(isat,:)];
JX=[JX ;jxe(isat,:)./rhoe(isat,:)];
JZI=[JZI ;jzi(isat,:)./rhoi(isat,:)];
JYI=[JYI ;jyi(isat,:)./rhoi(isat,:)];
JXI=[JXI ;jxi(isat,:)./rhoi(isat,:)];
NE=[NE ;rhoe(isat,:)];
NI=[NI ;rhoi(isat,:)];
end
SSX=EEY.*BBZ-EEZ.*BBY;
SSY=EEZ.*BBX-EEX.*BBZ;
SSZ=EEX.*BBY-EEY.*BBX;
Ndetrend=100;
SSXdet=SSX-tsmovavg(SSX,'s',Ndetrend);
BB=sqrt(BBX.^2+BBY.^2+BBZ.^2);
BBdet=BB-tsmovavg(BB,'s',Ndetrend);


figure(1)
pcolor((1:n1)*Dt*B0x,zzp,EEX)
shading interp
colorbar
title(['EX    x= ' num2str(xp(isat(1)))  '   y=   ' num2str(yp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('z')
set(gcf,'Renderer','zbuffer');
print('-dpng','EEXsat_xy.png' )
figure(2)
pcolor((1:n1)*Dt*B0x,zzp,EEY)
shading interp
colorbar
title(['EY    x= ' num2str(xp(isat(1)))  '   y=   ' num2str(yp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('z')
set(gcf,'Renderer','zbuffer');
print('-dpng','EEYsat_xy.png' )
figure(3)
pcolor((1:n1)*Dt*B0x,zzp,EEZ)
shading interp
colorbar
title(['EZ    x= ' num2str(xp(isat(1)))  '   y=   ' num2str(yp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('z')
set(gcf,'Renderer','zbuffer');
print('-dpng','EEZsat_xy.png' )

figure(4)
pcolor((1:n1)*Dt*B0x,zzp,BBX)
shading interp
colorbar
title(['BX    x= ' num2str(xp(isat(1)))  '   y=   ' num2str(yp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('z')
set(gcf,'Renderer','zbuffer');
print('-dpng','BBXsat_xy.png' )
figure(5)
pcolor((1:n1)*Dt*B0x,zzp,BBY)
shading interp
colorbar
title(['BY    x= ' num2str(xp(isat(1)))  '   y=   ' num2str(yp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('z')
set(gcf,'Renderer','zbuffer');
print('-dpng','BBYsat_xy.png' )
figure(6)
pcolor((1:n1)*Dt*B0x,zzp,BBZ)
shading interp
colorbar
title(['BZ    x= ' num2str(xp(isat(1)))  '   y=   ' num2str(yp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('z')
set(gcf,'Renderer','zbuffer');
print('-dpng','BBZsat_xy.png' )

figure(7)
pcolor((1:n1)*Dt*B0x,zzp,NE)
shading interp
colorbar
title(['\rho_e   x= ' num2str(xp(isat(1)))  '   y=   ' num2str(yp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('z')
set(gcf,'Renderer','zbuffer');
print('-dpng','NEsat_xy.png' )

figure(8)
pcolor((1:n1)*Dt*B0x,zzp,JPAR)
shading interp
colorbar
title(['u_{e||}   x= ' num2str(xp(isat(1)))  '   y=   ' num2str(yp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('z')
set(gcf,'Renderer','zbuffer');
print('-dpng','UEPARsat_xy.png' )

figure(9)
pcolor((1:n1)*Dt*B0x,zzp,JX)
shading interp
colorbar
title(['u_{ex}   x= ' num2str(xp(isat(1)))  '   y=   ' num2str(yp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('z')
set(gcf,'Renderer','zbuffer');
print('-dpng','UEXsat_xy.png' )
close(9)

figure(9)
pcolor((1:n1)*Dt*B0x,zzp,JY)
shading interp
colorbar
title(['u_{ex}   x= ' num2str(xp(isat(1)))  '   y=   ' num2str(yp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('z')
set(gcf,'Renderer','zbuffer');
print('-dpng','UEYsat_xy.png' )
close(9)

figure(10)
pcolor((1:n1)*Dt*B0x,zzp,JZ)
shading interp
colorbar
title(['u_{ez}   x= ' num2str(xp(isat(1)))  '   y=   ' num2str(yp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('z')
set(gcf,'Renderer','zbuffer');
print('-dpng','UEZsat_xy.png' )

figure(9) 
pcolor((1:n1)*Dt*B0x,zzp,JXI)
shading interp 
colorbar 
title(['u_{ix}   x= ' num2str(xp(isat(1)))  '   y=   ' num2str(yp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('z') 
set(gcf,'Renderer','zbuffer');
print('-dpng','UIXsat_xy.png' )
close(9) 
 
figure(9) 
pcolor((1:n1)*Dt*B0x,zzp,JYI)
shading interp 
colorbar 
title(['u_{ix}   x= ' num2str(xp(isat(1)))  '   y=   ' num2str(yp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('z') 
set(gcf,'Renderer','zbuffer');
print('-dpng','UIYsat_xy.png' )
close(9) 

figure(9) 
pcolor((1:n1)*Dt*B0x,zzp,JZI)
shading interp 
colorbar 
title(['u_{ix}   x= ' num2str(xp(isat(1)))  '   y=   ' num2str(yp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('z') 
set(gcf,'Renderer','zbuffer');
print('-dpng','UIZsat_xy.png' )
close(9) 

figure(11)
pcolor((1:n1)*Dt*B0x,zzp,NI)
shading interp
colorbar
title(['\rho_i   x= ' num2str(xp(isat(1)))  '   y=   ' num2str(yp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('z')
set(gcf,'Renderer','zbuffer');
print('-dpng','NIsat_xy.png' )

figure(12)
pcolor((1:n1)*Dt*B0x,zzp,SSX)
shading interp
colorbar
title(['ExB_x   x= ' num2str(xp(isat(1)))  '   y=   ' num2str(yp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('z')
set(gcf,'Renderer','zbuffer');
print('-dpng','SXsat_xy.png' )
close(12)
pcolor((1:n1)*Dt*B0x,zzp,SSXdet)
shading interp
colorbar
title(['ExB_x_{det}   x= ' num2str(xp(isat(1)))  '   y=   ' num2str(yp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('z')
set(gcf,'Renderer','zbuffer');
print('-dpng','SXsat_detrended_xy.png' )


figure(13)
pcolor((1:n1)*Dt*B0x,zzp,SSY)
shading interp
colorbar
title(['ExB_y   x= ' num2str(xp(isat(1)))  '   y=   ' num2str(yp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('z')
set(gcf,'Renderer','zbuffer');
print('-dpng','SYsat_xy.png' )

figure(14)
pcolor((1:n1)*Dt*B0x,zzp,SSZ)
shading interp
colorbar
title(['ExB_z   x= ' num2str(xp(isat(1)))  '   y=   ' num2str(yp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('z')
set(gcf,'Renderer','zbuffer');
print('-dpng','SZsat_xy.png' )

figure(15)
for i=1:1
y=SSX(i,:);
Fs = 1/Dt;                    % Sampling frequency
T = 1/Fs;                     % Sample time
L = max(size(y));             % Length of signal
t = (0:L-1)*T;                % Time vector
%y=sin(2*pi*t);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
hold off
loglog(f*2*pi,2*abs(Y(1:round(NFFT/2+1))),'k')
xlabel('\omega=2\pi f')
hold on
loglog([wpi wpi],[1e-10, 1e-5],'g')
loglog(sqrt(mratio)*[wpi wpi],[1e-10, 1e-5],'g--')
%loglog(sqrt(mratio)*[wpi_cavity wpi_cavity],[1e-6, 1e0],'g--')
loglog([wlh wlh],[1e-10, 1e-5],'r')
loglog([wce wce],[1e-10, 1e-5],'b--')
loglog([wci wci],[1e-10, 1e-5],'b')
loglog([1/Dt 1/Dt]*pi,[1e-10, 1e-5],'m')
%legend('spectrum','\omega_{pi}','\omega_{pe}','\omega_{pe,cavity}','\omega_{lh}','\omega_{ce}','\omega_{ci}')
legend('spectrum','\omega_{pi}','\omega_{pe}','\omega_{lh}','\omega_{ce}','\omega_{ci}','\pi/Dt')
title(['x= ' num2str(xp(isat(1)))  '   y=   ' num2str(yp(isat(1))) '   z=   ' num2str(zp(isat(1)))])
end
set(gcf,'Renderer','zbuffer');
print('-dpng', '-r300','SPECsat_xy.png' )

close(15)
figure(15)
y=SSXdet(i,Ndetrend+1:end);
Fs = 1/Dt;                    % Sampling frequency
T = 1/Fs;                     % Sample time
L = max(size(y));             % Length of signal
t = (0:L-1)*T;                % Time vector
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
hold off
loglog(f*2*pi,2*abs(Y(1:round(NFFT/2+1))),'k')
xlabel('\omega=2\pi f')
hold on
loglog([wpi wpi],[1e-14, 1e-8],'g')
loglog(sqrt(mratio)*[wpi wpi],[1e-14, 1e-8],'g--')
%loglog(sqrt(mratio)*[wpi_cavity wpi_cavity],[1e-6, 1e0],'g--')
loglog([wlh wlh],[1e-14, 1e-8],'r')
loglog([wce wce],[1e-14, 1e-8],'b--')
loglog([wci wci],[1e-14, 1e-8],'b')
loglog([1/Ndetrend/Dt 1/Ndetrend/Dt]*2*pi,[1e-14, 1e-8],'m')
%legend('spectrum','\omega_{pi}','\omega_{pe}','\omega_{pe,cavity}','\omega_{lh}','\omega_{ce}','\omega_{ci}')
legend('spectrum','\omega_{pi}','\omega_{pe}','\omega_{lh}','\omega_{ce}','\omega_{ci}','\omega_{detrend}','location','SouthWest')
title(['x= ' num2str(xp(isat(1)))  '   y=   ' num2str(yp(isat(1))) '   z=   ' num2str(zp(isat(1)))])

close(15)
figure(15)
y=BBdet(i,Ndetrend+1:end);
Fs = 1/Dt;                    % Sampling frequency
T = 1/Fs;                     % Sample time
L = max(size(y));             % Length of signal
t = (0:L-1)*T;                % Time vector
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
hold off
loglog(f*2*pi,2*abs(Y(1:round(NFFT/2+1))),'k')
xlabel('\omega=2\pi f')
hold on
loglog([wpi wpi],[1e-14, 1e-8],'g')
loglog(sqrt(mratio)*[wpi wpi],[1e-14, 1e-8],'g--')
%loglog(sqrt(mratio)*[wpi_cavity wpi_cavity],[1e-6, 1e0],'g--')
loglog([wlh wlh],[1e-14, 1e-8],'r')
loglog([wce wce],[1e-14, 1e-8],'b--')
loglog([wci wci],[1e-14, 1e-8],'b')
loglog([1/Ndetrend/Dt 1/Ndetrend/Dt]*2*pi,[1e-14, 1e-8],'m')
%legend('spectrum','\omega_{pi}','\omega_{pe}','\omega_{pe,cavity}','\omega_{lh}','\omega_{ce}','\omega_{ci}')
legend('spectrum','\omega_{pi}','\omega_{pe}','\omega_{lh}','\omega_{ce}','\omega_{ci}','\omega_{detrend}','location','SouthWest')
title(['x= ' num2str(xp(isat(1)))  '   y=   ' num2str(yp(isat(1))) '   z=   ' num2str(zp(isat(1)))])

set(gcf,'Renderer','zbuffer');
print('-dpng', '-r300','SPECBBsat_detrended_xy.png' )


figure(8)
[sx,sy] = meshgrid(15,1:1:10);   
ut=ones(size(JPAR));
h=streamline(stream2((1:n1)*Dt*B0x,zzp,ut*Dt*B0x,JZ,sx,sy,[1]));
set(h,'Color','white');
[sx,sy] = meshgrid(13,1:1:10);   
ut=ones(size(JPAR));
h=streamline(stream2((1:n1)*Dt*B0x,zzp,ut*Dt*B0x,JZ,sx,sy,[1]));
set(h,'Color','m');
set(gcf,'Renderer','zbuffer');
print('-dpng','Flowtraces_xy.png' )
