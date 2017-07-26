addpath 'matlab-parsek'

close all

nojump=1
if(nojump) 
vthr=.045;

%results_dir='/shared/gianni/laila13/'
results_dir='./'

filename=[results_dir 'settings.hdf']
Lx=hdf5read(filename,'/collective/Lx')
Ly=hdf5read(filename,'/collective/Ly')
B0x=hdf5read(filename,'/collective/B0x')
%B0x=hdf5read(filename,'/collective/Bx0')
Dt=hdf5read(filename,'/collective/Dt')
Dx=hdf5read(filename,'/collective/Dx')
XLEN=hdf5read(filename,'/topology/XLEN')
YLEN=hdf5read(filename,'/topology/YLEN')
mratio=abs(hdf5read(filename,'/collective/species_0/qom'))
Nprocs=hdf5read(filename,'/topology/Nprocs')

%ipx=XLEN/2-2; ipy=YLEN/2+2;
% for run 99 use: %ipy=YLEN/2-3;
ipy=round(YLEN/2);

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

n1=floor(n/16)

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
jxi=reshape(jxi(1:n1*16),16,n1);
jyi=reshape(jyi(1:n1*16),16,n1);
jzi=reshape(jzi(1:n1*16),16,n1);
edotj=(ex.*(jxe+jxi)+ey.*(jye+jyi)+ez.*(jze+jzi));
t=linspace(0,n1,n1);


%n1sat=40737;
n1sat=20000;
%n1sat=35000;

%isatx=1:2;isaty=1;
isatx=1:2;isaty=ISATY;
isat=(isatx-1)*8+isaty;
xxp=[xxp xp(isat)];
yyp=[yyp yp(isat)];
EEX=[EEX ;ex(isat,1:n1sat)];
EEY=[EEY ;ey(isat,1:n1sat)];
EEZ=[EEZ ;ez(isat,1:n1sat)];
BBX=[BBX ;bx(isat,1:n1sat)];
BBY=[BBY ;by(isat,1:n1sat)];
BBZ=[BBZ ;bz(isat,1:n1sat)];
JX=[JX ;jxe(isat,1:n1sat)./rhoe(isat,1:n1sat)];
JY=[JY ;jye(isat,1:n1sat)./rhoe(isat,1:n1sat)];
JZ=[JZ ;jze(isat,1:n1sat)./rhoe(isat,1:n1sat)];
NE=[NE ;rhoe(isat,1:n1sat)];
EDOTJ=[EDOTJ;edotj(isat,1:n1sat)];
end

npmax=max(size(EEX));

Dxs=xxp(2)-xxp(1);
AZ=-cumsum(BBY).*Dxs;
AZmin=min(AZ);
AZmax=max(AZ);
recon=AZmax-AZmin;

savva = [(1:npmax)'*Dt*B0x recon'/B0x];
save -ascii recon_sat.txt savva;   

plot((1:npmax)*Dt*B0x,recon(1:npmax)/B0x)
print -dpng recon-flux
recrate=diff(recon)/Dt;
plot((1.5:npmax-.5)*Dt*B0x,recrate(1:npmax-1)/B0x^2)
print -dpng recon-rate

recrateavg=tsmovavg(recrate,'s',100);
plot((1.5:npmax-.5)*Dt*B0x,recrateavg(1:npmax-1)/B0x^2)
print -dpng recon-rate-avg

SSX=EEY.*BBZ-EEZ.*BBY;
SSY=EEZ.*BBX-EEX.*BBZ;
SSZ=EEX.*BBY-EEY.*BBX;

%compute ohm non ideal terms (including hte hall term)
OHMX = EEX + JY.*BBZ - JZ.*BBY;
OHMY = EEY - JX.*BBZ + JZ.*BBX;
OHMZ = EEZ + JX.*BBY - JY.*BBX;


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

n1=n1sat;
callplot_sat_yz
