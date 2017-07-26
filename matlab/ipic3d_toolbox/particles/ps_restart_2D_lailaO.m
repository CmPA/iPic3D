clear all
addpath '/home/gianni/matlab2/matlab-parsek'
addpath '/home/gianni/matlab2'

results_dir='/shared02/gianni/lailaO1/data4/';
Nspecies=4
results_dir='/shared02/gianni/jean4/data4/';
Nspecies=1
results_dir='/shared02/gianni/shaped3/data2/';
Nspecies=4

filename=[results_dir 'settings.hdf'];
XLEN=double(hdf5read(filename,'/topology/XLEN'))
YLEN=double(hdf5read(filename,'/topology/YLEN'))
%ZLEN=1
ZLEN=double(hdf5read(filename,'/topology/ZLEN'))
Lx=double(hdf5read(filename,'/collective/Lx'))
Ly=double(hdf5read(filename,'/collective/Ly'))
%Lz=1
Lz=double(hdf5read(filename,'/collective/Lz'))
%B0x=double(hdf5read(filename,'/collective/B0x'))
B0x=double(hdf5read(filename,'/collective/Bx0'))
Nxc=double(hdf5read(filename,'/collective/Nxc'))
Nyc=double(hdf5read(filename,'/collective/Nyc'))
%Nzc=1
Nzc=double(hdf5read(filename,'/collective/Nzc'))
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

xminone=1e10;
xmaxone=-1e10;

ipx=round(XLEN/2);
% For LailaO
%for ipy=8:11
for ipy=YLEN/2
binnone=[];
% For LailaO
%for ipx=0*XLEN/4:XLEN/2
for ipx=0:XLEN
%ipy=8
%ipy=10

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
[xg,yg]=ndgrid(xmin:dx:xmax,ymin:dy:ymax);

%close all

vmax_graph=[.25 .06 .25 .06]/3;

%Run Laila0
%is=3
%is2=3
is=2
is2=is
% For LailaO
%vnorm=.3
vnorm=.1/sqrt(256);
% for shaped3
vnorm=.1
is=2
is2=is

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
Bx=squeeze(Bx(1,:,:))';
By=squeeze(By(1,:,:))';
Bz=squeeze(Bz(1,:,:))';


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

second_species=0;
if(second_species)
q=[q ;squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is2).Groups(1).Datasets(it)))];
x=[x ;squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is2).Groups(5).Datasets(it)))];
y=[y ;squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is2).Groups(6).Datasets(it)))];
%z=[z ;squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is2).Groups(7).Datasets(it)))];
u=[u ;squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is2).Groups(2).Datasets(it)))];
v=[v ;squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is2).Groups(3).Datasets(it)))];
w=[w ;squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is2).Groups(4).Datasets(it)))];
end

ii=y<mean(y);
x=x(ii);
y=y(ii);
q=q(ii);
u=u(ii);
v=v(ii);
w=w(ii);

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

[totnum,nbinu,xrange,urange]=spaziofasi(x',u',q',0,-vnorm,vnorm);
%[totnum,nbinu,xrange,urange]=spaziofasi(uperp1',upar',q',0,-vnorm,vnorm);
x1=xmin;x2=xmax;
imagesc([x1 x2],urange',log10(nbinu(1:end-1,:))')
binnone=[binnone; nbinu(1:end-1,:)];
xminone=min(xminone,min(x));
xmaxone=max(xmaxone,max(x));
%imagesc([x1 x2],urange',(nbinu(1:end-1,:))')
%set(gca,'xdir','reverse','TickDir','out')
xlabel('x/d_i','fontsize',14)
ylabel('Vpar/c','fontsize',14)
title(['y=' num2str(min(y)) ' , ' num2str(max(y))])
set(gca,'fontsize',[14])
set(gcf,'Renderer','zbuffer');
print('-dpng','-r300',['SP' num2str(is) 'phasespace2stream' num2str(ipx) '.png'])
%print('-depsc','-painters',['SP' num2str(is) 'phasespace2stream' num2str(ipx) '.eps'])
figure(2)
 bar(urange',mean(nbinu)./sum(mean(nbinu)))
title(['Proc(iy)=' num2str(ipx) '  x=[' num2str(min(x)) ' , ' num2str(max(x)) ']  y=[' num2str(min(y)) ' , ' num2str(max(y)),']'])
ylabel('f(Vpar)','fontsize',14)
xlabel('Vpar/c','fontsize',14)
set(gca,'fontsize',[14])
set(gcf,'Renderer','zbuffer');
['SP' num2str(is) 'distr2stream' num2str(ipx) '.png']
print('-dpng','-r300',['SP' num2str(is) 'distr2stream' num2str(ipx) '.png'])
print('-depsc','-painters',['SP' num2str(is) 'distr2stream' num2str(ipx) '.eps'])


end

%pause(.1)
end

scrsz = get(0,'ScreenSize')
close all
figure('Position',[1 scrsz(4)/2 scrsz(3)*4 scrsz(4)/2])
imagesc([xminone xmaxone],[-vnorm vnorm],(log10(binnone')))
hold on
plot([xminone xmaxone],[0 0],'w','linewidth',[2])
set(gca,'ydir','reverse','TickDir','out') 
axis xy
xlabel('x/d_i','fontsize',18)
ylabel('Vpar/c','fontsize',18)
title(['y=' num2str(min(y)) ' , ' num2str(max(y))],'fontsize',18)
set(gca,'fontsize',[18])
set(gcf,'Renderer','zbuffer');

saveas(gcf,['aaa' num2str(ipy) '.fig'])
set(gcf,'PaperUnits', 'inches')
xwidth=24; ywidth=6;
set(gcf,'PaperPosition', [0 0 xwidth ywidth])
saveas(gcf,['aaa' num2str(ipy) '.png'],'png')
%print('-dpng','-r300',['aaa' num2str(ipy) '.png'])
end
