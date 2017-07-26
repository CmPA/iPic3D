
close all
clear all

results_dir='/shared/gianni/tred54/'

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

ipx=XLEN/2-3%-4; 
ipy=YLEN/2+1; %+1
%ipy=YLEN/2+2*0;

ipz=ZLEN/2;
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
satxp=mean(xp)
satyp=mean(yp)
satzp=mean(zp)


a=fscanf(fid,'%f',[14 inf])';

fclose(fid);
skip=0;
ez=a(:,6+skip);
ey=a(:,5+skip);
ex=a(:,4+skip);
bz=a(:,3+skip);
by=a(:,2+skip);
bx=a(:,1+skip);
rhoe=-a(:,13+skip)*4*pi;
rhoi=a(:,14+skip)*4*pi;

[n m]=size(ez);
n1=floor(n/27);

ex=reshape(ex(1:n1*27),27,n1);
ey=reshape(ey(1:n1*27),27,n1);
ez=reshape(ez(1:n1*27),27,n1);
bx=reshape(bx(1:n1*27),27,n1);
by=reshape(by(1:n1*27),27,n1);
bz=reshape(bz(1:n1*27),27,n1);
rhoe=reshape(rhoe(1:n1*27),27,n1);
rhoi=reshape(rhoi(1:n1*27),27,n1);

isat=1 %7
ex=ex(isat,:);
ey=ey(isat,:);
ez=ez(isat,:);
bx=bx(isat,:);
by=by(isat,:);
bz=bz(isat,:);
rhoe=rhoe(isat,:);
rhoi=rhoi(isat,:);

Navg=50/2

bxavg=tsmovavg(bx,'s',Navg);
byavg=tsmovavg(by,'s',Navg);
bzavg=tsmovavg(bz,'s',Navg);
bavg=sqrt(bxavg.^2+byavg.^2+bzavg.^2)

dbx=bx-bxavg;
dby=by-byavg;
dbz=bz-bzavg;

dex=ex-tsmovavg(ex,'s',Navg);
dey=ey-tsmovavg(ey,'s',Navg);
dez=ez-tsmovavg(ez,'s',Navg);

depar=(dex.*bxavg + dey.*byavg + dez.*bzavg)./bavg;

bmag2D=sqrt(bxavg.^2+byavg.^2);
deperp1=(byavg.*dex-bxavg.*dey)./bmag2D;

perp2x=bzavg.*bxavg./(bavg.*bmag2D);
perp2y=bzavg.*byavg./(bavg.*bmag2D);
perp2z=-bmag2D./bavg;
deperp2=perp2x.*dex+perp2y.*dey+perp2z.*dez;


dbpar=(dbx.*bxavg + dby.*byavg + dbz.*bzavg)./bavg;

bmag2D=sqrt(bxavg.^2+byavg.^2);
dbperp1=(byavg.*dbx-bxavg.*dby)./bmag2D;

perp2x=bzavg.*bxavg./(bavg.*bmag2D);
perp2y=bzavg.*byavg./(bavg.*bmag2D);
perp2z=-bmag2D./bavg;
dbperp2=perp2x.*dbx+perp2y.*dby+perp2z.*dbz;


ex=dex;
ey=dey;
ez=dez;
Bmean=mean(sqrt(bx.^2+by.^2+bz.^2))
%time=B0x*Dt*(1:max(size(bx)));
time=Bmean*Dt*(1:max(size(bx)));

h=figure(1)
set(h,'Position', [608 41 1118 962])
subplot(4,1,1)
plot(time, rhoe, time, rhoi,time, rhoi-rhoe)
xlabel('\omega_{ci}t')
ylabel('rhoe')
legend('rhoe','rhoi','rhonet','location','EastOutside')
title(['x=' num2str(xp(isat)-Lx/2) '   y=' num2str(yp(isat)-Ly/2) '   z=' num2str(zp(isat)-Lz/2) ])
subplot(4,1,2)
plot(time, bx,time,by,time,bz)
xlabel('\omega_{ci}t')
ylabel('B')
legend('Bx','By','Bz','location','EastOutside')
subplot(4,1,3)
plot(time, dbx,time,dby,time,dbz)
xlabel('\omega_{ci}t')
ylabel('dB')
legend('dBx','dBy','dBz','location','EastOutside')
subplot(4,1,4)
plot(time, ex,time,ey,time,ez)
xlabel('\omega_{ci}t')
ylabel('E')
legend('Ex','Ey','Ez','location','EastOutside')


figure(2)
subplot(2,2,1)
plot(ex,ey)
xlabel('Ex')
ylabel('Ey')
subplot(2,2,2)
plot(ex,ez)
xlabel('Ex')
ylabel('Ez')
subplot(2,2,3)
plot(ex,ey)
xlabel('Ey')
ylabel('Ez')

figure(3)
i=sqrt(-1);
subplot(3,1,1)
teta=angle(ex+i*ey)
plot(time,teta)
xlabel('\omega_{ci}t')
ylabel('Angle(Ex,Ey)')
subplot(3,1,2)
teta=angle(ex+i*ez)
plot(time,teta)
xlabel('\omega_{ci}t')
ylabel('Angle(Ex,Ez)')
subplot(3,1,3)
teta=angle(ey+i*ez)
plot(time,teta)
xlabel('\omega_{ci}t')
ylabel('Angle(Ey,Ez)')


figure(4)
i=sqrt(-1);
subplot(3,1,1)
plot(time,ex./sqrt(ex.^2+ey.^2))
xlabel('\omega_{ci}t')
ylabel('Ex/Exy')
ylim([-1 1])
subplot(3,1,2)
plot(time,ex./sqrt(ex.^2+ez.^2))
xlabel('\omega_{ci}t')
ylabel('Ex/Exz')
ylim([-1 1])
subplot(3,1,3)
plot(time,ey./sqrt(ey.^2+ez.^2))
xlabel('\omega_{ci}t')
ylabel('Ey/(Ezy)')
ylim([-2 2])

figure(5)
i=sqrt(-1);
subplot(3,1,1)
teta=angle(dbx+i*dby)
plot(time,teta)
xlabel('\omega_{ci}t')
ylabel('Angle(dBx,dBy)')
subplot(3,1,2)
teta=angle(dbx+i*dbz)
plot(time,teta)
xlabel('\omega_{ci}t')
ylabel('Angle(dBx,dBz)')
subplot(3,1,3)
teta=angle(dby+i*dbz)
plot(time,teta)
xlabel('\omega_{ci}t')
ylabel('Angle(dBy,dBz)')


figure(6)
i=sqrt(-1);
subplot(3,1,1)
plot(time,dbx./sqrt(dbx.^2+dby.^2))
xlabel('\omega_{ci}t')
ylabel('dBx/dBxy')
ylim([-1 1])
subplot(3,1,2)
plot(time,dbx./sqrt(dbx.^2+dbz.^2))
xlabel('\omega_{ci}t')
ylabel('dBx/dBxz')
ylim([-1 1])
subplot(3,1,3)
plot(time,dby./sqrt(dby.^2+dbz.^2))
xlabel('\omega_{ci}t')
ylabel('dBy/(dBzy)')
ylim([-2 2])

figure(7)
ii=time>2.9 & time<3.15;
subplot(3,2,1)
plot(time(ii),depar(ii))
xlabel('\omega_{ci}t')
ylabel('E_{||}')
subplot(3,2,3)
plot(time(ii),dbperp1(ii))
xlabel('\omega_{ci}t')
ylabel('B_{\perp1}')
subplot(3,2,5)
plot(time(ii),dbperp2(ii))
xlabel('\omega_{ci}t')
ylabel('B_{\perp2}')
subplot(3,2,[2 4 6])
plot(dbperp1(ii),dbperp2(ii))
axis equal
xlabel('E_{\perp 1}')
ylabel('E_{\perp 2}')
print -dpng hodogram
print -depsc -painters hodogram
figure(8)
plot3(time(ii),dbperp1(ii),dbperp2(ii))
xlabel('\omega_{ci}t')
ylabel('B_{\perp 1}')
zlabel('B_{\perp 2}')
title(['xsat=' num2str(xp(1)) ', ysat=',num2str(yp(1)) ', zsat=',num2str(zp(1)) ])
print -dpng hodogram3D
print -depsc -painters hodogram3D
