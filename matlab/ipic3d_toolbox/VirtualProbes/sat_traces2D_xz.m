addpath 'matlab-parsek'

close all

nojump=1
if(nojump) 
vthr=.045;

results_dir='./';

filename=[results_dir 'settings.hdf']
Lx=hdf5read(filename,'/collective/Lx')
Ly=hdf5read(filename,'/collective/Ly')
B0x=hdf5read(filename,'/collective/B0x')
Dt=hdf5read(filename,'/collective/Dt')
XLEN=hdf5read(filename,'/topology/XLEN')
YLEN=hdf5read(filename,'/topology/YLEN')
mratio=abs(hdf5read(filename,'/collective/species_0/qom'))
Nprocs=hdf5read(filename,'/topology/Nprocs')


ipx=XLEN/2+15;
%ipx=2;


yyp=[];
EEX=[];
EEY=[];
EEZ=[];
BBX=[];
BBY=[];
BBZ=[];
JPAR=[];
NE=[];
for ipy=1:YLEN;
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
%n1=14000;
%n1=50000;

ex=reshape(ex(1:n1*16),16,n1);
ey=reshape(ey(1:n1*16),16,n1);
ez=reshape(ez(1:n1*16),16,n1);
bx=reshape(bx(1:n1*16),16,n1);
by=reshape(by(1:n1*16),16,n1);
bz=reshape(bz(1:n1*16),16,n1);
rhoe=reshape(rhoe(1:n1*16),16,n1);
jepar=reshape(jepar(1:n1*16),16,n1);
t=linspace(0,n1,n1);


isatx=1;isaty=1:8;
isat=(isatx-1)*8+isaty;
yyp=[yyp yp(isat)];
EEX=[EEX ;ex(isat,:)];
EEY=[EEY ;ey(isat,:)];
EEZ=[EEZ ;ez(isat,:)];
BBX=[BBX ;bx(isat,:)];
BBY=[BBY ;by(isat,:)];
BBZ=[BBZ ;bz(isat,:)];
JPAR=[JPAR ;jepar(isat,:)./rhoe(isat,:)];
NE=[NE ;rhoe(isat,:)];
end


SSX=EEY.*BBZ-EEZ.*BBY;
SSY=EEZ.*BBX-EEX.*BBZ;
SSZ=EEX.*BBY-EEY.*BBX;



npmax=max(size(EEX));
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


n1=npmax;

callplot_sat_xz

return

figure(1)
pcolor((1:n1)*Dt*B0x,yyp,EEX)
shading interp
colorbar
title(['EX   x=   ' num2str(xp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('y')
set(gcf,'Renderer','zbuffer');
print('-dpng','EEXsat_xz.png' )
figure(2)
pcolor((1:n1)*Dt*B0x,yyp,EEY)
shading interp
colorbar
title(['EY    x=   ' num2str(xp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('y')
set(gcf,'Renderer','zbuffer');
print('-dpng','EEYsat_xz.png' )
figure(3)
pcolor((1:n1)*Dt*B0x,yyp,EEZ)
shading interp
colorbar
title(['EZ    x=   ' num2str(xp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('y')
set(gcf,'Renderer','zbuffer');
print('-dpng','EEZsat_xz.png' )

figure(4)
pcolor((1:n1)*Dt*B0x,yyp,BBX)
shading interp
colorbar
title(['BX    x=   ' num2str(xp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('y')
set(gcf,'Renderer','zbuffer');
print('-dpng','BBXsat_xz.png' )
figure(5)
pcolor((1:n1)*Dt*B0x,yyp,BBY)
shading interp
colorbar
title(['BY    x=   ' num2str(xp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('y')
set(gcf,'Renderer','zbuffer');
print('-dpng','BBYsat_xz.png' )
figure(6)
pcolor((1:n1)*Dt*B0x,yyp,BBZ)
shading interp
colorbar
title(['BZ    x=   ' num2str(xp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('y')
set(gcf,'Renderer','zbuffer');
print('-dpng','BBZsat_xz.png' )

figure(7)
pcolor((1:n1)*Dt*B0x,yyp,NE)
shading interp
colorbar
title(['\rho_e   x=   ' num2str(xp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('y')
set(gcf,'Renderer','zbuffer');
print('-dpng','NEsat_xz.png' )

figure(8)
pcolor((1:n1)*Dt*B0x,yyp,JPAR)
shading interp
colorbar
title(['u_e   x=   ' num2str(xp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('y')
set(gcf,'Renderer','zbuffer');
print('-dpng','UEsat_xz.png' )

figure(12)
pcolor((1:n1)*Dt*B0x,yyp,SSX)
shading interp
colorbar
title(['ExB_x      x=   ' num2str(xp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('y')
set(gcf,'Renderer','zbuffer');
print('-dpng','SXsat_xz.png' )

figure(13)
pcolor((1:n1)*Dt*B0x,yyp,SSY)
shading interp
colorbar
title(['ExB_y     x=   ' num2str(xp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('y')
set(gcf,'Renderer','zbuffer');
print('-dpng','SYsat_xz.png' )

figure(14)
pcolor((1:n1)*Dt*B0x,yyp,SSZ)
shading interp
colorbar
title(['ExB_z     x=   ' num2str(xp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('y')
set(gcf,'Renderer','zbuffer');
print('-dpng','SZsat_xz.png' )

figure(15)
pcolor((Ndetrend+1:n1)*Dt*B0x,yyp,SSXdet(:,Ndetrend+1:n1))
shading interp
colorbar
title(['[Detrended(E)xDetrendend(B)]_x      x=   ' num2str(xp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('y')
set(gcf,'Renderer','zbuffer');
print('-dpng','SXsatdet_xz.png' )

figure(16)
pcolor((2*Ndetrend+1:n1)*Dt*B0x,yyp,SSXdetavg(:,2*Ndetrend+1:n1))
shading interp
colorbar
title(['<[Detrended(E)xDetrendend(B)]>_x      x=   ' num2str(xp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('y')
set(gcf,'Renderer','zbuffer');
print('-dpng','SXsatdetavg_xz.png' )


