addpath 'matlab-parsek'

close all

vthr=.045;

results_dir='/shared/gianni/tred43/';

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
%ipz=ZLEN/2; ipy=YLEN/2-1;
ipz=ZLEN/2; ipy=YLEN/2;

xxp=[];
EEX=[];
EEY=[];
EEZ=[];
BBX=[];
BBY=[];
BBZ=[];
UPAR=[];
UX=[];
UZ=[];
NE=[];
for ipx=1:XLEN;

ip=(ipx-1)*YLEN*ZLEN+(ipy-1)*ZLEN+ipz-1;

nome=[results_dir 'VirtualSatelliteTraces' num2str(ip) '.txt']

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
%n1=14555


%n1=14000;

ex=reshape(ex(1:n1*27),27,n1);
ey=reshape(ey(1:n1*27),27,n1);
ez=reshape(ez(1:n1*27),27,n1);
bx=reshape(bx(1:n1*27),27,n1);
by=reshape(by(1:n1*27),27,n1);
bz=reshape(bz(1:n1*27),27,n1);
rhoe=reshape(rhoe(1:n1*27),27,n1);
jepar=reshape(jepar(1:n1*27),27,n1);
jxe=reshape(jxe(1:n1*27),27,n1);
jze=reshape(jze(1:n1*27),27,n1);
t=linspace(0,n1,n1);



isatx=1:3;isaty=1;isatz=1;
isat=(isatx-1)*3*3+(isaty-1)*3+isatz;
%[xp(isat);yp(isat);zp(isat)]
xxp=[xxp xp(isat)]
EEX=[EEX ;ex(isat,:)];
EEY=[EEY ;ey(isat,:)];
EEZ=[EEZ ;ez(isat,:)];
BBX=[BBX ;bx(isat,:)];
BBY=[BBY ;by(isat,:)];
BBZ=[BBZ ;bz(isat,:)];
UPAR=[UPAR ;jepar(isat,:)./rhoe(isat,:)];
UX=[UX ;jxe(isat,:)./rhoe(isat,:)];
UZ=[UZ ;jze(isat,:)./rhoe(isat,:)];
NE=[NE ;rhoe(isat,:)];
end

n1p=10000; %time before restart second portion

results_dir='/shared/gianni/tred43.2/';
xxp=[];
EEX2=[];
EEY2=[];
EEZ2=[];
BBX2=[];
BBY2=[];
BBZ2=[];
UPAR2=[];
UX2=[];
UZ2=[];
NE2=[];
for ipx=1:XLEN;

ip=(ipx-1)*YLEN*ZLEN+(ipy-1)*ZLEN+ipz-1;

nome=[results_dir 'VirtualSatelliteTraces' num2str(ip) '.txt']

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
%n1=14555


%n1=14000;

ex=reshape(ex(1:n1*27),27,n1);
ey=reshape(ey(1:n1*27),27,n1);
ez=reshape(ez(1:n1*27),27,n1);
bx=reshape(bx(1:n1*27),27,n1);
by=reshape(by(1:n1*27),27,n1);
bz=reshape(bz(1:n1*27),27,n1);
rhoe=reshape(rhoe(1:n1*27),27,n1);
jepar=reshape(jepar(1:n1*27),27,n1);
jxe=reshape(jxe(1:n1*27),27,n1);
jze=reshape(jze(1:n1*27),27,n1);

t=linspace(0,n1,n1);



isatx=1:3;isaty=1;isatz=1;
isat=(isatx-1)*3*3+(isaty-1)*3+isatz;
%[xp(isat);yp(isat);zp(isat)]
xxp=[xxp xp(isat)]
EEX2=[EEX2 ;ex(isat,:)];
EEY2=[EEY2 ;ey(isat,:)];
EEZ2=[EEZ2 ;ez(isat,:)];
BBX2=[BBX2 ;bx(isat,:)];
BBY2=[BBY2 ;by(isat,:)];
BBZ2=[BBZ2 ;bz(isat,:)];
UPAR2=[UPAR2 ;jepar(isat,:)./rhoe(isat,:)];
UX2=[UX2 ;jxe(isat,:)./rhoe(isat,:)];
UZ2=[UZ2 ;jze(isat,:)./rhoe(isat,:)];
NE2=[NE2 ;rhoe(isat,:)];
end

n1=n1+n1p;

EEX=[EEX(:,1:n1p) EEX2];
EEY=[EEY(:,1:n1p) EEY2];
EEZ=[EEZ(:,1:n1p) EEZ2];
BBX=[BBX(:,1:n1p) BBX2];
BBY=[BBY(:,1:n1p) BBY2];
BBZ=[BBZ(:,1:n1p) BBZ2];
UPAR=[UPAR(:,1:n1p) UPAR2];
UX=[UX(:,1:n1p) UX2];
UZ=[UZ(:,1:n1p) UZ2];
NE=[NE(:,1:n1p) NE2];

npmax=15000;

EEX=EEX(:,1:npmax);
EEY=EEY(:,1:npmax);
EEZ=EEZ(:,1:npmax);
BBX=BBX(:,1:npmax);
BBY=BBY(:,1:npmax);
BBZ=BBZ(:,1:npmax);
UPAR=UPAR(:,1:npmax);
UX=UX(:,1:npmax);
UZ=UZ(:,1:npmax);
NE=NE(:,1:npmax);
	
	
SSX=EEY.*BBZ-EEZ.*BBY;
SSY=EEZ.*BBX-EEX.*BBZ;
SSZ=EEX.*BBY-EEY.*BBX;

figure(1)
pcolor((1:npmax)*Dt*B0x,xxp,EEX)
shading interp
colorbar
title(['EX    z= ' num2str(zp(isatz))  '   y=   ' num2str(yp(isaty))])
xlabel('\omega_{ci}t')
ylabel('x')
set(gcf,'Renderer','zbuffer');
print('-dpng','EEXsat_yz.png' )
saveas(gcf,'EEXsat_yz.fig')

figure(2)
pcolor((1:npmax)*Dt*B0x,xxp,EEY)
shading interp
colorbar
title(['EY    z= ' num2str(zp(isatz))  '   y=   ' num2str(yp(isaty))])
xlabel('\omega_{ci}t')
ylabel('x')
set(gcf,'Renderer','zbuffer');
print('-dpng','EEYsat_yz.png' )
saveas(gcf,'EEYsat_yz.fig')

figure(3)
pcolor((1:npmax)*Dt*B0x,xxp,EEZ)
shading interp
colorbar
title(['EZ    z= ' num2str(zp(isatz))  '   y=   ' num2str(yp(isaty))])
xlabel('\omega_{ci}t')
ylabel('x')
set(gcf,'Renderer','zbuffer');
print('-dpng','EEZsat_yz.png' )
saveas(gcf,'EEZsat_yz.fig')


figure(4)
pcolor((1:npmax)*Dt*B0x,xxp,BBX)
shading interp
colorbar
title(['BX    z= ' num2str(zp(isatz))  '   y=   ' num2str(yp(isaty))])
xlabel('\omega_{ci}t')
ylabel('x')
set(gcf,'Renderer','zbuffer');
print('-dpng','BBXsat_yz.png' )
saveas(gcf,'BBXsat_yz.fig')

figure(5)
pcolor((1:npmax)*Dt*B0x,xxp,BBY)
shading interp
colorbar
title(['BY    z= ' num2str(zp(isatz))  '   y=   ' num2str(yp(isaty))])
xlabel('\omega_{ci}t')
ylabel('x')
set(gcf,'Renderer','zbuffer');
print('-dpng','BBYsat_yz.png' )
saveas(gcf,'BBYsat_yz.fig')

figure(6)
pcolor((1:npmax)*Dt*B0x,xxp,BBZ)
shading interp
colorbar
title(['BZ    z= ' num2str(zp(isatz))  '   y=   ' num2str(yp(isaty))])
xlabel('\omega_{ci}t')
ylabel('x')
set(gcf,'Renderer','zbuffer');
print('-dpng','BBZsat_yz.png' )
saveas(gcf,'BBZsat_yz.fig')

figure(7)
pcolor((1:npmax)*Dt*B0x,xxp,NE)
shading interp
colorbar
title(['\rho_e   z= ' num2str(zp(isatz))  '   y=   ' num2str(yp(isaty))])
xlabel('\omega_{ci}t')
ylabel('x')
set(gcf,'Renderer','zbuffer');
print('-dpng','NEsat_yz.png' )
saveas(gcf,'NEsat_yz.fig')

figure(8)
pcolor((1:npmax)*Dt*B0x,xxp,UPAR)
shading interp
colorbar
title(['u_e   z= ' num2str(zp(isatz))  '   y=   ' num2str(yp(isaty))])
xlabel('\omega_{ci}t')
ylabel('x')
set(gcf,'Renderer','zbuffer');
print('-dpng','UPARsat_yz.png' )
saveas(gcf,'UPARsat_yz.fig')

figure(9)
pcolor((1:npmax)*Dt*B0x,xxp,UX)
shading interp
colorbar
title(['u_{xe}   z= ' num2str(zp(isatz))  '   y=   ' num2str(yp(isaty))])
xlabel('\omega_{ci}t')
ylabel('x')
set(gcf,'Renderer','zbuffer');
print('-dpng','UXsat_yz.png' )
saveas(gcf,'UXsat_yz.fig')

figure(10)
pcolor((1:npmax)*Dt*B0x,xxp,UZ)
shading interp
colorbar
title(['u_{ze}   z= ' num2str(zp(isatz))  '   y=   ' num2str(yp(isaty))])
xlabel('\omega_{ci}t')
ylabel('x')
set(gcf,'Renderer','zbuffer');
print('-dpng','UYsat_yz.png' )
saveas(gcf,'UYsat_yz.fig')


figure(12)
pcolor((1:npmax)*Dt*B0x,xxp,SSX)
shading interp
colorbar
title(['ExB_x   z= ' num2str(zp(isat(1)))  '   y=   ' num2str(yp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('x')
set(gcf,'Renderer','zbuffer');
print('-dpng','SXsat_yz.png' )
saveas(gcf,'SXsat_yz.fig')

figure(13)
pcolor((1:npmax)*Dt*B0x,xxp,SSY)
shading interp
colorbar
title(['ExB_y   z= ' num2str(zp(isat(1)))  '   y=   ' num2str(yp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('z')
set(gcf,'Renderer','zbuffer');
print('-dpng','SYsat_yz.png' )
saveas(gcf,'SYsat_yz.fig')

figure(14)
pcolor((1:npmax)*Dt*B0x,xxp,SSZ)
shading interp
colorbar
title(['ExB_z   z= ' num2str(zp(isat(1)))  '   y=   ' num2str(yp(isat(1)))])
xlabel('\omega_{ci}t')
ylabel('z')
set(gcf,'Renderer','zbuffer');
print('-dpng','SZsat_yz.png' )
saveas(gcf,'SZsat_yz.fig')
