addpath 'matlab-parsek'

close all

vthr=.045;

results_dir='/shared/gianni/tred48/';
results_dir2='/shared/gianni/tred48.2/';

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

ipx=XLEN/2+3;
ipy=YLEN/2-1;
ipz=ZLEN/2;
ip=(ipx-1)*YLEN*ZLEN+(ipy-1)*ZLEN+ipz;

nome=[results_dir 'VirtualSatelliteTraces' num2str(ip) '.txt']

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

[n m]=size(bx);

n0=mean(rhoi-rhoe)/2;
b0=sqrt(mean(bx.^2+by.^2+bz.^2));
wci=b0;
wpi=1*sqrt(n0);
wce=wci*mratio;
wlh=1/sqrt(1/wce/wci+1/wpi^2);
%wpi=1 %apparently the plasma oscillations are generated elsewhere where n0=1

n1=floor(n/27)
%n1=14000;
bx=reshape(bx(1:n1*27),27,n1);
by=reshape(by(1:n1*27),27,n1);
bz=reshape(bz(1:n1*27),27,n1);
ex=reshape(ex(1:n1*27),27,n1);
ey=reshape(ey(1:n1*27),27,n1);
ez=reshape(ez(1:n1*27),27,n1);
rhoe=reshape(rhoe(1:n1*27),27,n1);
t=linspace(0,n1,n1);
b=sqrt(bx.*bx+by.*by+bz.*bz);
epar=(ex.*bx+ey.*by+ez.*bz)./b;




figure(1)
plot(t,ex(1,:),'r')
title(['xsat=' num2str(xp(1)) '  ysat=',num2str(yp(1))])
ylabel('Ex')
xlabel('\omega_{pi}t')
hold on 


nome=[results_dir2 'VirtualSatelliteTraces' num2str(ip) '.txt']

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

[n m]=size(bx);

n0=mean(rhoi-rhoe)/2;
b0=sqrt(mean(bx.^2+by.^2+bz.^2));
wci=b0;
wpi=1*sqrt(n0);
wce=wci*mratio;
wlh=1/sqrt(1/wce/wci+1/wpi^2);
%wpi=1 %apparently the plasma oscillations are generated elsewhere where n0=1

n2=floor(n/27)
nrest=12000;
bx=reshape(bx(1:n2*27),27,n2);
by=reshape(by(1:n2*27),27,n2);
bz=reshape(bz(1:n2*27),27,n2);
ex=reshape(ex(1:n2*27),27,n2);
ey=reshape(ey(1:n2*27),27,n2);
ez=reshape(ez(1:n2*27),27,n2);
rhoe=reshape(rhoe(1:n2*27),27,n2);
t=linspace(nrest,nrest+n2,n2);
b=sqrt(bx.*bx+by.*by+bz.*bz);
epar=(ex.*bx+ey.*by+ez.*bz)./b;


figure(1)
plot(t,ex(1,:),'b')
