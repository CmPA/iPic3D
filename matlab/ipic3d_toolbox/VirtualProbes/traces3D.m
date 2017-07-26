addpath 'matlab-parsek'

close all
clear all


results_dir='/shared/gianni/tred54.2/';
%results_dir='/shared/gianni/run100.2/';
results_dir='/data1/gianni/HRmaha3D3/traces/';


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
vthr=hdf5read(filename,'/collective/species_0/uth')
Nprocs=hdf5read(filename,'/topology/Nprocs')

ipx=XLEN/2-3;
ipy=YLEN/2+1;
ipz=ZLEN/2;
ip=(ipx-1)*YLEN*ZLEN+(ipy-1)*ZLEN+ipz;

nome=[results_dir 'VirtualSatelliteTraces' num2str(ip) '.txt']
%nome=[results_dir 'VirtualSatelliteTraces1876.txt']

if(exist(nome)~=2)
disp(['transfer ' nome])
return
end

fid=fopen(nome);
for i=1:27
x=fscanf(fid,'%f',3); 
xp(i)=x(1);
yp(i)=x(2);
zp(i)=x(3);
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

n0=mean(rhoi-rhoe)/2;
b0=sqrt(mean(bx.^2+by.^2+bz.^2));
wci0=B0x;
wci=b0;
wpi=1*sqrt(n0);
wce=wci*mratio;
wlh=1/sqrt(1/wce/wci+1/wpi^2);
%wpi=1 %apparently the plasma oscillations are generated elsewhere where n0=1

[n m]=size(bx);
n2=floor(n/27)
bx=reshape(bx(1:n2*27),27,n2);
by=reshape(by(1:n2*27),27,n2);
bz=reshape(bz(1:n2*27),27,n2);
ex=reshape(ex(1:n2*27),27,n2);
ey=reshape(ey(1:n2*27),27,n2);
ez=reshape(ez(1:n2*27),27,n2);
rhoe=reshape(rhoe(1:n2*27),27,n2);
b=sqrt(bx.*bx+by.*by+bz.*bz);
epar=(ex.*bx+ey.*by+ez.*bz)./b;

Navg=20;

[exdet,dummy2]=analisi_traces(ex,b,'Ex',Dt,Navg,wpi,wci0,wci,wlh,mratio,xp,yp,zp);
[eydet,dummy2]=analisi_traces(ey,b,'Ey',Dt,Navg,wpi,wci0,wci,wlh,mratio,xp,yp,zp);
[ezdet,dummy2]=analisi_traces(ez,b,'Ez',Dt,Navg,wpi,wci0,wci,wlh,mratio,xp,yp,zp);
[dummy1,dummy2]=analisi_traces(bx,b,'Bx',Dt,Navg,wpi,wci0,wci,wlh,mratio,xp,yp,zp);
[dummy1,dummy2]=analisi_traces(by,b,'By',Dt,Navg,wpi,wci0,wci,wlh,mratio,xp,yp,zp);
[dummy1,dummy2]=analisi_traces(bz,b,'Bz',Dt,Navg,wpi,wci0,wci,wlh,mratio,xp,yp,zp);
[dummy1,dummy2]=analisi_traces(rhoe,b,'Rhoe',Dt,Navg,wpi,wci0,wci,wlh,mratio,xp,yp,zp);


%dby=diff(by,1,2)/Dt;
%dby=dby(:,Navg-1:end);

faisvd=0
if(faisvd)
[u,s,v]=svd(by');
end
%plot(u(:,1)*s(1,1)*max(v(1,:)))
%title('topos')
if(faisvd)
plot(wci0*t,byavg(1,:), wci0*t,u(:,1)*s(1,1)*max(v(1,:)),'r')
title(['xsat=' num2str(xp(1)) '  ysat=',num2str(yp(1))])
end
