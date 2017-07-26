clear all
addpath '/home/gianni/matlab2/matlab-parsek'
addpath '/home/gianni/matlab2'

results_dir='/shared02/gianni/maha2/data4/';
Nspecies=2
species_plot=1

filename=[results_dir 'settings.hdf'];
XLEN=double(hdf5read(filename,'/topology/XLEN'))
YLEN=double(hdf5read(filename,'/topology/YLEN'))
ZLEN=1
%ZLEN=double(hdf5read(filename,'/topology/ZLEN'))
Lx=double(hdf5read(filename,'/collective/Lx'))
Ly=double(hdf5read(filename,'/collective/Ly'))
Lz=1
%Lz=double(hdf5read(filename,'/collective/Lz'))
%B0x=double(hdf5read(filename,'/collective/B0x'))
B0x=double(hdf5read(filename,'/collective/Bx0'))
Nxc=double(hdf5read(filename,'/collective/Nxc'))
Nyc=double(hdf5read(filename,'/collective/Nyc'))
Nzc=1
%Nzc=double(hdf5read(filename,'/collective/Nzc'))
vth=double(hdf5read(filename,'/collective/species_0/uth'))

% use for pqrsek 2D 
%B0x=double(hdf5read(filename,'/collective/B0x'))

Dt=double(hdf5read(filename,'/collective/Dt'))

mp = 1.6726e-27;
me= mp/256;
e= 1.6022e-19;
mu0=4*pi*1e-7;
eps0=8.8542*1.e-12;
c=1/sqrt(mu0*eps0);
TeoTi=1/5;
energy_code=me*c^2/e/TeoTi;


%ipx=1;
%for  ipy=[2 11]
for ipx=1:XLEN
%ipy=11
ipy=YLEN/2

if(ZLEN>1)
ipz=ZLEN/2;
else
ipz=1;
end

ip=(ipx-1)*YLEN*ZLEN+(ipy-1)*ZLEN+ipz-1;

xmin=double(0+(ipx-1)*Lx/XLEN)
xmax=double(xmin+Lx/XLEN)
ymin=double(0+(ipy-1)*Ly/YLEN)
ymax=double(ymin+double(Ly/YLEN))
zmin=double(0+(ipz-1)*Lz/ZLEN)
zmax=double(zmin+double(Lz/ZLEN))

nx=Nxc/XLEN+1;
ny=Nyc/YLEN+1;
nz=Nzc/ZLEN+1;

dx=double(xmax-xmin)/double(nx-1);
dy=double(ymax-ymin)/double(ny-1);
dz=double(zmax-zmin)/double(nz-1);

[yg,xg]=ndgrid(ymin:dy:ymax,xmin:dx:xmax);


%close all

vmax_graph=[.25 .06 .25 .06]/3;

is=2

vmax = vmax_graph(is);
partfile='part';
partfile='restart';
if exist([results_dir partfile num2str(ip) '.hdf']) ==2 

info=hdf5info([results_dir partfile num2str(ip) '.hdf']);
[nt dummy]=size(info.GroupHierarchy.Groups(3).Groups(is).Groups(5).Datasets(:));


%infof=hdf5info([results_dir 'proc' num2str(ip) '.hdf']);

%[ntf dummy]=size(infof.GroupHierarchy.Groups(2).Groups(4).Datasets(:));

it=nt
%ex=squeeze(hdf5read(infof.GroupHierarchy.Groups(2).Groups(4).Datasets(it)));   
%ex=squeeze(hdf5read(infof.GroupHierarchy.Groups(1).Groups(4).Datasets(it)));   


Bx=squeeze(hdf5read(info.GroupHierarchy.Groups(1).Groups(1).Datasets(it)));
By=squeeze(hdf5read(info.GroupHierarchy.Groups(1).Groups(2).Datasets(it)));
Bz=squeeze(hdf5read(info.GroupHierarchy.Groups(1).Groups(3).Datasets(it)));
% for runs with ipic3D
Bx=squeeze(Bx(1,:,:));
By=squeeze(By(1,:,:));
Bz=squeeze(Bz(1,:,:));


%
% 1 = q
% 2 = u
% 3 = v
% 4 = w
% 5 = x
% 6 = y
% 7 = z

q=squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is).Groups(1).Datasets(it)));
x=squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is).Groups(5).Datasets(it)));
y=squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is).Groups(6).Datasets(it)));
%z=squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is).Groups(7).Datasets(it)));
u=squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is).Groups(2).Datasets(it)));
v=squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is).Groups(3).Datasets(it)));
w=squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is).Groups(4).Datasets(it)));
cutout=0
if(cutout)
ii=y>mean(y(:));
%ii=y>0;
q=q(ii);
x=x(ii);
y=y(ii);
%z=z(ii);
u=u(ii);
v=v(ii);
w=w(ii);
ii=y<mean(y(:));
%ii=y>0;
q=q(ii);
x=x(ii);
y=y(ii);
%z=z(ii);
u=u(ii);
v=v(ii);
w=w(ii);
end
bxp=interpn(xg,yg,Bx,x,y);
byp=interpn(xg,yg,By,x,y);
bzp=interpn(xg,yg,Bz,x,y);
bp=sqrt(bxp.^2+byp.^2+bzp.^2);
bp2D=sqrt(bxp.^2+byp.^2);
perp2x=bzp.*bxp./(bp.*bp2D);
perp2y=bzp.*byp./(bp.*bp2D);
perp2z=-bp2D./bp;


upar=(u.*bxp+v.*byp+w.*bzp)./(bp+1e-10);
uperp1=(byp.*u-bxp.*v)./bp2D;

uperp2=perp2x.*u+perp2y.*v+perp2z.*w;

epar=.5*upar.^2;
eperp=.5*(uperp1.^2+uperp2.^2);

epar=epar*energy_code;
eperp=eperp*energy_code;


figure(1)
[totnum,nbinu,xrange,urange]=spaziofasi(x',upar',q',0);
x1=-min(xrange(:))/Lx*30-15;
x2=-max(xrange(:))/Lx*30-15;
imagesc([x1 x2],urange',log10(nbinu(1:end-1,:))')
set(gca,'xdir','reverse','TickDir','out')
xlabel('x/R_E','fontsize',14)
ylabel('Vpar/c','fontsize',14)
set(gca,'fontsize',[14])
set(gcf,'Renderer','zbuffer');
print('-dpng','-r300',['SP' num2str(is) 'phasespace2stream' num2str(ipx) '.png'])
print('-depsc','-painters',['SP' num2str(is) 'phasespace2stream' num2str(ipx) '.eps'])
figure(2)
 bar(urange',mean(nbinu)./sum(mean(nbinu)))
title(['Proc(iy)=' num2str(ipy) '  x(GSM)=[' num2str(-min(x/Lx)*30-15) ' , ' num2str(-max(x/Lx)*30-15) ']  z(GSM)=[' num2str(min(y/Ly)*12-9) ' , ' num2str(max(y/Ly)*12-9),']'])
ylabel('f(Vpar)','fontsize',14)
xlabel('Vpar/c','fontsize',14)
set(gca,'fontsize',[14])
set(gcf,'Renderer','zbuffer');
print('-dpng','-r300',['SP' num2str(is) 'distr2stream' num2str(ipx) '.png'])
print('-depsc','-painters',['SP' num2str(is) 'distr2stream' num2str(ipx) '.eps'])

end

pause(.1)
end
