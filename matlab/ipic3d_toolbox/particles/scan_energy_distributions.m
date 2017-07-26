addpath '/home/gianni/matlab2/matlab-parsek'
addpath '/home/gianni/matlab2'

% This works i the run has saved the particles
%results_dir='/shared/gianni/lhdi6redux/';
results_dir='/shared02/gianni/tred60/data1/';
%results_dir='/shared/gianni/laila11b/';

filename=[results_dir 'settings.hdf'];
XLEN=double(hdf5read(filename,'/topology/XLEN'))
YLEN=double(hdf5read(filename,'/topology/YLEN'))
ZLEN=double(hdf5read(filename,'/topology/ZLEN'))
Lx=double(hdf5read(filename,'/collective/Lx'))
Ly=double(hdf5read(filename,'/collective/Ly'))
Lz=double(hdf5read(filename,'/collective/Lz'))
B0x=double(hdf5read(filename,'/collective/Bx0'))
Nxc=double(hdf5read(filename,'/collective/Nxc'))
Nyc=double(hdf5read(filename,'/collective/Nyc'))
Nzc=double(hdf5read(filename,'/collective/Nzc'))
vth=double(hdf5read(filename,'/collective/species_0/uth'))

% use for pqrsek 2D 
%B0x=double(hdf5read(filename,'/collective/B0x'))

Dt=double(hdf5read(filename,'/collective/Dt'))

bins=0:(vth/5).^2:(5*vth)^2;
histpar=[];
histperp=[];
xhist=[];
bhist=[];
bxhist=[];
byhist=[];
bzhist=[];


is=1 % for electrons
% is=2 % for ions

for  ipx=1:XLEN

ipy=YLEN/2+2;
ipz=ZLEN/2;

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

[xg,yg,zg]=ndgrid(xmin:dx:xmax,ymin:dy:ymax,zmin:dz:zmax);

close all

vmax_graph=[.25 .06 .25 .06];

vmax = vmax_graph(is);

partfile='restart';

info=hdf5info([results_dir partfile num2str(ip) '.hdf']);
[it dummy]=size(info.GroupHierarchy.Groups(3).Groups(is).Groups(5).Datasets(:));


infof=hdf5info([results_dir 'proc' num2str(ip) '.hdf']);

ex=squeeze(hdf5read(infof.GroupHierarchy.Groups(1).Groups(4).Datasets(it)));   
Bx=squeeze(hdf5read(info.GroupHierarchy.Groups(1).Groups(1).Datasets(it)));
By=squeeze(hdf5read(info.GroupHierarchy.Groups(1).Groups(2).Datasets(it)));
Bz=squeeze(hdf5read(info.GroupHierarchy.Groups(1).Groups(3).Datasets(it)));


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
z=squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is).Groups(7).Datasets(it)));

bxp=interpn(xg,yg,zg,permute(Bx,[3 2 1]),x,y,z);
byp=interpn(xg,yg,zg,permute(By,[3 2 1]),x,y,z);
bzp=interpn(xg,yg,zg,permute(Bz,[3 2 1]),x,y,z);
bp=sqrt(bxp.^2+byp.^2+bzp.^2);
bp2D=sqrt(bxp.^2+byp.^2);
perp2x=bzp.*bxp./(bp.*bp2D);
perp2y=bzp.*byp./(bp.*bp2D);
perp2z=-bp2D./bp;

u=squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is).Groups(2).Datasets(it)));
v=squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is).Groups(3).Datasets(it)));
w=squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is).Groups(4).Datasets(it)));

upar=(u.*bxp+v.*byp+w.*bzp)./(bp+1e-10);
uperp1=(byp.*u-bxp.*v)./bp2D;

uperp2=perp2x.*u+perp2y.*v+perp2z.*w;

epar=.5*upar.^2;
eperp=.5*(uperp1.^2+uperp2.^2);


histpar=[histpar; hist(epar,bins)];
histperp=[histperp; hist(eperp,bins)]; 
xhist=[xhist;mean(x)];
bhist=[bhist;mean(bp)];
bxhist=[bxhist;mean(bxp)];
byhist=[byhist;mean(byp)];
bzhist=[bzhist;mean(bzp)];

end


figure
subplot(2,1,1)
for i=1:16
	histperp2(i,:) = histperp(i,:).*bins;
end					 
pcolor(xhist,bins,log(histperp2'))
xlabel('x/d_i')
ylabel('energy')
title('Perp Energy Flux Ef(E)')
shading interp
subplot(2,1,2)
plot(xhist,bhist,xhist,bxhist,xhist,byhist,xhist,bzhist)
				legend('B','B_x','B_y','B_z')
set(gcf,'Renderer','zbuffer');
print('-dpng','EnergyFluxPerp.png')


figure
subplot(2,1,1)
for i=1:16
	histpar2(i,:) = histpar(i,:).*bins;
end					 
pcolor(xhist,bins,log(histpar2'))
xlabel('x/d_i')
ylabel('energy')
title('Par Energy Flux Ef(E)')
shading interp
subplot(2,1,2)
plot(xhist,bhist,xhist,bxhist,xhist,byhist,xhist,bzhist)
					  legend('B','B_x','B_y','B_z')
set(gcf,'Renderer','zbuffer');
print('-dpng','EnergyFluxPar.png')

