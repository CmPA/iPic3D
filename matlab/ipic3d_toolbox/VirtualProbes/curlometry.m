addpath 'matlab-parsek'

close all
clear all

vthr=.045;

results_dir='';

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

ipx = XLEN/2 -2
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
for ipz=1:1;
ip=(ipx-1)*YLEN*ZLEN+(ipy-1)*ZLEN+ipz-1;

ip=913
nome=[results_dir 'VirtualSatelliteTraces' num2str(ip) '.txt']
%system(['gunzip ' nome])
if(exist(nome)~=2)
system(['scp giovanni@128.97.76.22:/rd/giovanni/HRmaha3D1/probes/data1/' 'VirtualSatelliteTraces' num2str(ip) '.txt' ' ' results_dir])
end
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

n1=4200

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

