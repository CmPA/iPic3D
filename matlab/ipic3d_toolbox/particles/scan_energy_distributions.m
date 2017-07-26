clear all 
close all

addpath(genpath('../../ipic3d_toolbox'));

% This works i the run has saved the particles
results_dir='/Users/gianni/Dropbox/Science/codes/iPic3D-github/build/data/';

global Lx Ly Lz Xgsmrange Ygsmrange Zgsmrange dx dy dz XLEN YLEN ZLEN initial_time Nx Ny Nz


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
qom=double(hdf5read(filename,'/collective/species_0/qom'))
% use for pqrsek 2D 
%B0x=double(hdf5read(filename,'/collective/B0x'))

Dt=double(hdf5read(filename,'/collective/Dt'))


Xgsmrange= [-39.8 -8];
Zgsmrange= [-10.8 5];
Ygsmrange= [-4 15.8];


is=2 % for electrons
% is=2 % for ions

if(is==1) 
    % set electron max speed
    vmax=.2
else
     % set ion max speed
    vmax=.2/sqrt(-qom);
end

bins=0:(vmax/25).^2:(vmax)^2;

histpar=[];
histperp=[];
xhist=[];
bhist=[];
bxhist=[];
byhist=[];
bzhist=[];

nsubx=20;
nsuby=1;
nsubz=1;


for  ipx=1:XLEN

ipy=YLEN/2;
ipz=ZLEN;

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

partfile='restart';

info=hdf5info([results_dir partfile num2str(ip) '.hdf']);
[it dummy]=size(info.GroupHierarchy.Groups(3).Groups(is).Groups(5).Datasets(:));
   
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

dxsub=(xmax-xmin)/nsubx;
dysub=(ymax-ymin)/nsuby;
dzsub=(zmax-zmin)/nsubz;

ixsub=floor((x-xmin)/dxsub);
iysub=floor((y-ymin)/dysub);
izsub=floor((z-zmin)/dzsub);

suby=0;subz=0;
for subx=0:nsubx-1
ii=(ixsub==subx)&(iysub==suby)&(izsub==subz);


hist_epar=hist(epar(ii),bins);
hist_eperp=hist(eperp(ii),bins);
histpar=[histpar; hist_epar];
histperp=[histperp; hist_eperp]; 
xhist=[xhist;mean(x(ii))];
bhist=[bhist;mean(bp(ii))];
bxhist=[bxhist;mean(bxp(ii))];
byhist=[byhist;mean(byp(ii))];
bzhist=[bzhist;mean(bzp(ii))];

end


end

[nscan, nbin]=size(histperp);


RE=6371000;
box=1.6671*RE;
npart=32*32*14*14;
n0=.25e6;
c=300000e3;
dv=mean(diff(sqrt(bins)))*c;
dx=Lx/XLEN/nsubx*c;


mp = 1.6726e-27;
me= mp/256;
e= 1.6022e-19;
mu0=4*pi*1e-7;
eps0=8.8542*1.e-12;
c=1/sqrt(mu0*eps0);
k=1.38e-23;
energy_code=me*c^2/e;


code_E = 2060.21;
code_B = 6.87213e-06;
code_B=code_B *1e9; % to convert from Tesla to nT
code_J = 1.20082e-05;
code_J = code_J*1e9; % to convert to nA/m^2
code_V = 2.99792e+08;
code_V=code_V/1e3; %to convert to Km/s
code_T = 1.50326e-10;
code_n = 0.25;
e=1.6e-19;

%convert to keV
TeoTi=1/5;
code_T=code_T/e/1e3/TeoTi;


histperp=histperp/box^3/dv/dx/npart*box^3*n0;
histpar=histpar/box^3/dv/dx/npart*box^3*n0;
bins=bins*energy_code;

figure
subplot(2,1,1)
for i=1:nscan
	histperp_flux(i,:) = histperp(i,:).*bins;
end					 
pcolor(gsmx(xhist),bins,log(histperp_flux'))
set(gca,'Xdir','reverse')
xlabel('x/d_i')
ylabel('energy [eV]')
title('Perp Energy Flux Ef(E)')
shading interp
subplot(2,1,2)
plot(gsmx(xhist),bhist,gsmx(xhist),bxhist,gsmx(xhist),byhist,gsmx(xhist),bzhist)
set(gca,'Xdir','reverse')
				legend('B','B_x','B_y','B_z')
set(gcf,'Renderer','zbuffer');
print('-dpng','EnergyFluxPerp.png')


figure
subplot(2,1,1)
for i=1:nscan
	histpar_flux(i,:) = histpar(i,:).*bins;
end					 
pcolor(gsmx(xhist),bins,log(histpar_flux'))
set(gca,'Xdir','reverse')
xlabel('x/d_i')
ylabel('energy [eV]')
title('Par Energy Flux Ef(E)')
shading interp
subplot(2,1,2)
plot(gsmx(xhist),bhist,gsmx(xhist),bxhist,gsmx(xhist),byhist,gsmx(xhist),bzhist)
set(gca,'Xdir','reverse')
					  legend('B','B_x','B_y','B_z')
set(gcf,'Renderer','zbuffer');
print('-dpng','EnergyFluxPar.png')

