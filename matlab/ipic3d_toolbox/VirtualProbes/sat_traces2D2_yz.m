addpath 'matlab-parsek'

close all

vthr=.045;

rootdir=pwd;
results_dir=[rootdir '/'];

filename=[results_dir 'settings.hdf']
Lx=hdf5read(filename,'/collective/Lx')
Ly=hdf5read(filename,'/collective/Ly')
B0x=hdf5read(filename,'/collective/B0x')
Dt=hdf5read(filename,'/collective/Dt')
XLEN=hdf5read(filename,'/topology/XLEN')
YLEN=hdf5read(filename,'/topology/YLEN')
mratio=abs(hdf5read(filename,'/collective/species_0/qom'))
Nprocs=hdf5read(filename,'/topology/Nprocs')


readin=1

if(readin) 
%ipx=XLEN/2-2; ipy=YLEN/2+2;
% for run 99 use: %ipy=YLEN/2-3;
ipy=YLEN/2-1;

xxp=[];
yyp=[];
EEX=[];
EEY=[];
EEZ=[];
BBX=[];
BBY=[];
BBZ=[];
JX=[];
JY=[];
JZ=[];
JiX=[];
JiY=[];
JiZ=[];
NE=[];
EDOTJ=[];
for ipx=1:XLEN;
ip=(ipx-1)*YLEN+(ipy-1);

nome=[results_dir 'VirtualSatelliteTraces' num2str(ip) '.txt']
system(['gunzip ' nome])

fid=fopen(nome);
for i=1:16
x=fscanf(fid,'%f',2); 
xp(i)=x(1);
yp(i)=x(2);
end


a=fscanf(fid,'%f',[14 inf])';

fclose(fid);
skip=0;
bx=a(:,1+skip);
by=a(:,2+skip);
bz=a(:,3+skip);
ex=a(:,4+skip);
ey=a(:,5+skip);
ez=a(:,6+skip);
jxe=a(:,7+skip)*4*pi;
jye=a(:,8+skip)*4*pi;
jze=a(:,9+skip)*4*pi;
jxi=a(:,10+skip)*4*pi;
jyi=a(:,11+skip)*4*pi;
jzi=a(:,12+skip)*4*pi;
rhoe=a(:,13+skip)*4*pi;
rhoi=a(:,14+skip)*4*pi;
edotj=(ex.*(jxe+jxi)+ey.*(jye+jyi)+ez.*(jze+jzi));
b=sqrt(bx.*bx+by.*by+bz.*bz);
epar=(ex.*bx+ey.*by+ez.*bz)./b;
jepar=(jxe.*bx+jye.*by+jze.*bz)./b;

[n m]=size(bx);

n0=mean(rhoi-rhoe)/2;
b0=sqrt(mean(bx.^2+by.^2+bz.^2));
wci=b0;
wpi=1*sqrt(n0);
wce=wci*mratio;
wlh=1/sqrt(1/wce/wci+1/wpi^2);
%wpi=1 %apparently the plasma oscillations are generated elsewhere where n0=1

n1=floor(n/16)
%n1=14000;

ex=reshape(ex(1:n1*16),16,n1);
ey=reshape(ey(1:n1*16),16,n1);
ez=reshape(ez(1:n1*16),16,n1);
bx=reshape(bx(1:n1*16),16,n1);
by=reshape(by(1:n1*16),16,n1);
bz=reshape(bz(1:n1*16),16,n1);
rhoe=reshape(rhoe(1:n1*16),16,n1);
jxe=reshape(jxe(1:n1*16),16,n1);
jye=reshape(jye(1:n1*16),16,n1);
jze=reshape(jze(1:n1*16),16,n1);
rhoi=reshape(rhoi(1:n1*16),16,n1);
jxi=reshape(jxi(1:n1*16),16,n1);
jyi=reshape(jyi(1:n1*16),16,n1);
jzi=reshape(jzi(1:n1*16),16,n1);
edotj=reshape(edotj(1:n1*16),16,n1);
t=linspace(0,n1,n1);

n1p=40000; %time before restart second portion
%n1p=50000; %time before restart second portion

isatx=1:4;isaty=ISATY;
isat=(isatx-1)*4+isaty;
yp(isat)
xxp=[xxp xp(isat)];
yyp=[yyp yp(isat)];
EEX=[EEX ;ex(isat,1:n1p)];
EEY=[EEY ;ey(isat,1:n1p)];
EEZ=[EEZ ;ez(isat,1:n1p)];
BBX=[BBX ;bx(isat,1:n1p)];
BBY=[BBY ;by(isat,1:n1p)];
BBZ=[BBZ ;bz(isat,1:n1p)];
JX=[JX ;jxe(isat,1:n1p)./rhoe(isat,1:n1p)];
JY=[JY ;jye(isat,1:n1p)./rhoe(isat,1:n1p)];
JZ=[JZ ;jze(isat,1:n1p)./rhoe(isat,1:n1p)];
JiX=[JiX ;jxi(isat,1:n1p)./rhoi(isat,1:n1p)];
JiY=[JiY ;jyi(isat,1:n1p)./rhoi(isat,1:n1p)];
JiZ=[JiZ ;jzi(isat,1:n1p)./rhoi(isat,1:n1p)];
NE=[NE ;rhoe(isat,1:n1p)];
EDOTJ=[EDOTJ;edotj(isat,1:n1p)];

end


results_dir=[rootdir '.2/'];

EEX2=[];
EEY2=[];
EEZ2=[];
BBX2=[];
BBY2=[];
BBZ2=[];
JX2=[];
JY2=[];
JZ2=[];
JiX2=[];
JiY2=[];
JiZ2=[];
NE2=[];
EDOTJ2=[];
for ipx=1:XLEN;
ip=(ipx-1)*YLEN+(ipy-1);

nome=[results_dir 'VirtualSatelliteTraces' num2str(ip) '.txt']
system(['gunzip ' nome])

fid=fopen(nome);
for i=1:16
x=fscanf(fid,'%f',2); 
end


a=fscanf(fid,'%f',[14 inf])';

fclose(fid);
skip=0;
bx=a(:,1+skip);
by=a(:,2+skip);
bz=a(:,3+skip);
ex=a(:,4+skip);
ey=a(:,5+skip);
ez=a(:,6+skip);
jxe=a(:,7+skip)*4*pi;
jye=a(:,8+skip)*4*pi;
jze=a(:,9+skip)*4*pi;
jxi=a(:,10+skip)*4*pi;
jyi=a(:,11+skip)*4*pi;
jzi=a(:,12+skip)*4*pi;
rhoe=a(:,13+skip)*4*pi;
rhoi=a(:,14+skip)*4*pi;
edotj=(ex.*(jxe+jxi)+ey.*(jye+jyi)+ez.*(jze+jzi));
b=sqrt(bx.*bx+by.*by+bz.*bz);
epar=(ex.*bx+ey.*by+ez.*bz)./b;
jepar=(jxe.*bx+jye.*by+jze.*bz)./b;

[n m]=size(bx);

n0=mean(rhoi-rhoe)/2;
b0=sqrt(mean(bx.^2+by.^2+bz.^2));
wci=b0;
wpi=1*sqrt(n0);
wce=wci*mratio;
wlh=1/sqrt(1/wce/wci+1/wpi^2);
%wpi=1 %apparently the plasma oscillations are generated elsewhere where n0=1

n1=floor(n/16)
%n1=14000;

ex=reshape(ex(1:n1*16),16,n1);
ey=reshape(ey(1:n1*16),16,n1);
ez=reshape(ez(1:n1*16),16,n1);
bx=reshape(bx(1:n1*16),16,n1);
by=reshape(by(1:n1*16),16,n1);
bz=reshape(bz(1:n1*16),16,n1);
rhoe=reshape(rhoe(1:n1*16),16,n1);
jxe=reshape(jxe(1:n1*16),16,n1);
jye=reshape(jye(1:n1*16),16,n1);
jze=reshape(jze(1:n1*16),16,n1);
rhoi=reshape(rhoi(1:n1*16),16,n1);
jxi=reshape(jxi(1:n1*16),16,n1);
jyi=reshape(jyi(1:n1*16),16,n1);
jzi=reshape(jzi(1:n1*16),16,n1);
t=linspace(0,n1,n1);
edotj=reshape(edotj(1:n1*16),16,n1);

n2p=20000; %run157, run159
%n2p=28000; %run158
n2p=22400; %run157

isatx=1:4;isaty=ISATY;
isat=(isatx-1)*4+isaty;
EEX2=[EEX2 ;ex(isat,1:n2p)];
EEY2=[EEY2 ;ey(isat,1:n2p)];
EEZ2=[EEZ2 ;ez(isat,1:n2p)];
BBX2=[BBX2 ;bx(isat,1:n2p)];
BBY2=[BBY2 ;by(isat,1:n2p)];
BBZ2=[BBZ2 ;bz(isat,1:n2p)];
JX2=[JX2 ;jxe(isat,1:n2p)./rhoe(isat,1:n2p)];
JY2=[JY2 ;jye(isat,1:n2p)./rhoe(isat,1:n2p)];
JZ2=[JZ2 ;jze(isat,1:n2p)./rhoe(isat,1:n2p)];
JiX2=[JiX2 ;jxi(isat,1:n2p)./rhoi(isat,1:n2p)];
JiY2=[JiY2 ;jyi(isat,1:n2p)./rhoi(isat,1:n2p)];
JiZ2=[JiZ2 ;jzi(isat,1:n2p)./rhoi(isat,1:n2p)];
NE2=[NE2 ;rhoe(isat,1:n2p)];
EDOTJ2=[EDOTJ2 ;edotj(isat,1:n2p)];

end


n1=n2p+n1p;

EEX=[EEX(:,1:n1p) EEX2];
clear EEX2;
EEY=[EEY(:,1:n1p) EEY2];
clear EEY2;
EEZ=[EEZ(:,1:n1p) EEZ2];
clear EEZ2;
BBX=[BBX(:,1:n1p) BBX2];
clear BBX2;
BBY=[BBY(:,1:n1p) BBY2];
clear BBY2;
BBZ=[BBZ(:,1:n1p) BBZ2];
clear BBZ2;
JX=[JX(:,1:n1p) JX2];
clear JX2;
JY=[JY(:,1:n1p) JY2];
clear JY2;
JZ=[JZ(:,1:n1p) JZ2];
clear JZ2;
JiX=[JiX(:,1:n1p) JiX2];
clear JiX2;
JiY=[JiY(:,1:n1p) JiY2];
clear JiY2;
JiZ=[JiZ(:,1:n1p) JiZ2];
clear JiZ2;
NE=[NE(:,1:n1p) NE2];
clear NE2;
EDOTJ=[EDOTJ(:,1:n1p) EDOTJ2];
clear EDOTJ2;

SSX=EEY.*BBZ-EEZ.*BBY;
SSY=EEZ.*BBX-EEX.*BBZ;
SSZ=EEX.*BBY-EEY.*BBX;

%compute ohm non ideal terms (including hte hall term)
OHMX = EEX + JY.*BBZ - JZ.*BBY;
OHMY = EEY - JX.*BBZ + JZ.*BBX;
OHMZ = EEZ + JX.*BBY - JY.*BBX;
OHMiX = EEX + JiY.*BBZ - JiZ.*BBY;
OHMiY = EEY - JiX.*BBZ + JiZ.*BBX;
OHMiZ = EEZ + JiX.*BBY - JiY.*BBX;


npmax=max(size(EEX));
if(detrended)
Ndetrend=100;
EEXdet=EEX(:,1:npmax)-tsmovavg(EEX(:,1:npmax),'s',Ndetrend);
EEYdet=EEY(:,1:npmax)-tsmovavg(EEY(:,1:npmax),'s',Ndetrend);
EEZdet=EEZ(:,1:npmax)-tsmovavg(EEZ(:,1:npmax),'s',Ndetrend);
BBXdet=BBX(:,1:npmax)-tsmovavg(BBX(:,1:npmax),'s',Ndetrend);
BBYdet=BBY(:,1:npmax)-tsmovavg(BBY(:,1:npmax),'s',Ndetrend);
BBZdet=BBZ(:,1:npmax)-tsmovavg(BBZ(:,1:npmax),'s',Ndetrend);

SSXdet=EEYdet.*BBZdet-EEZdet.*BBYdet;

SSXdetavg=tsmovavg(SSXdet(:,1:npmax),'s',Ndetrend);
tsmovavg(EEX(:,1:npmax),'s',Ndetrend);
end

end

n1=npmax;

callplot_sat_yz



