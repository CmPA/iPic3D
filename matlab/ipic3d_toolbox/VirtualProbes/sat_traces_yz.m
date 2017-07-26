addpath 'matlab-parsek'

close all

vthr=.045;

results_dir='./'

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

%ipx=XLEN/2-2; ipy=YLEN/2+2;
ipz=ZLEN/2; ipy=YLEN/2-1;
%tred54
ipz = ZLEN*2/10; ipy=YLEN*10/15;%9/15;%8/15;%was 11/15

%whistler in DF
ipy = YLEN/2 +1
ipz = ZLEN/2


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
mean(xp);
mean(yp);
mean(zp);


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

n0=mean(rhoi-rhoe)/2;
b0=sqrt(mean(bx.^2+by.^2+bz.^2));
wci=b0;
wpi=1*sqrt(n0);
wce=wci*mratio;
wlh=1/sqrt(1/wce/wci+1/wpi^2);
%wpi=1 %apparently the plasma oscillations are generated elsewhere where n0=1


n1=floor(n/27);


ex=reshape(ex(1:n1*27),27,n1);
ey=reshape(ey(1:n1*27),27,n1);
ez=reshape(ez(1:n1*27),27,n1);
bx=reshape(bx(1:n1*27),27,n1);
by=reshape(by(1:n1*27),27,n1);
bz=reshape(bz(1:n1*27),27,n1);
rhoe=reshape(rhoe(1:n1*27),27,n1);
jxe=reshape(jxe(1:n1*27),27,n1);
jye=reshape(jye(1:n1*27),27,n1);
jze=reshape(jze(1:n1*27),27,n1);
jxi=reshape(jxi(1:n1*27),27,n1);
jyi=reshape(jyi(1:n1*27),27,n1);
jzi=reshape(jzi(1:n1*27),27,n1);
edotj=(ex.*(jxe+jxi)+ey.*(jye+jyi)+ez.*(jze+jzi));
t=linspace(0,n1,n1);

%n1sat=12000;
n1sat=3200

isatx=1:3;isaty=ISATY;isatz=1;
isat=(isatx-1)*3*3+(isaty-1)*3+isatz;
%[xp(isat);yp(isat);zp(isat)]
xxp=[xxp xp(isat)]
yyp=[yyp yp(isat)]
EEX=[EEX ;ex(isat,1:n1sat)];
EEY=[EEY ;ey(isat,1:n1sat)];
EEZ=[EEZ ;ez(isat,1:n1sat)];
BBX=[BBX ;bx(isat,1:n1sat)];
BBY=[BBY ;by(isat,1:n1sat)];
BBZ=[BBZ ;bz(isat,1:n1sat)];
JX=[JX ;jxe(isat,1:n1sat)./rhoe(isat,1:n1sat)];
JY=[JY ;jye(isat,1:n1sat)./rhoe(isat,1:n1sat)];
JZ=[JZ ;jze(isat,1:n1sat)./rhoe(isat,1:n1sat)];
JiX=[JiX ;jxi(isat,1:n1sat)];
JiY=[JiY ;jyi(isat,1:n1sat)];
JiZ=[JiZ ;jzi(isat,1:n1sat)];
NE=[NE ;rhoe(isat,1:n1sat)];
EDOTJ=[EDOTJ;edotj(isat,1:n1sat)];
end
SSX=EEY.*BBZ-EEZ.*BBY;
SSY=EEZ.*BBX-EEX.*BBZ;
SSZ=EEX.*BBY-EEY.*BBX;

%compute ohm non ideal terms (including hte hall term)
OHMX = EEX + JY.*BBZ - JZ.*BBY;
OHMY = EEY - JX.*BBZ + JZ.*BBX;
OHMZ = EEZ + JX.*BBY - JY.*BBX;

npmax=max(size(EEX));
n1=n1sat;
callplot_sat_yz
