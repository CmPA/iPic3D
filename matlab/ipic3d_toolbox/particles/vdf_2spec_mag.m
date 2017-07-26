addpath(genpath('../../ipic3d_toolbox'));

global Lx Ly Lz Xgsmrange Ygsmrange Zgsmrange dx dy dz XLEN YLEN ZLEN initial_time Nx Ny Nz
close all

results_dir='/shared02/gianni/tred60/data1/';
%results_dir='/nobackup/glapenta/1mar08C/data3/';
results_dir='/Users/gianni/Documents/iPic3D-open/data/';

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
vth_ele=double(hdf5read(filename,'/collective/species_0/uth'))
vth_ion=double(hdf5read(filename,'/collective/species_1/uth'))


Dt=double(hdf5read(filename,'/collective/Dt'))


Xgsmrange= [-39.8 -8];
Zgsmrange= [-10.8 5];
Ygsmrange= [-4 15.8];


ipx=XLEN/2;
for  ipx=1:XLEN

ipy=YLEN/2+1;
ipz=ZLEN/2+1;
ipz=1

%1 mar08C
%ipy=7;ipz=10;
%Ygsm=11
%ipz = round(ZLEN * (Ygsm-Ygsmrange(1)) / (Ygsmrange(2)-Ygsmrange(1)) )
%Zgsm=1
%ipy = round(YLEN * (Zgsm-Zgsmrange(1)) / (Zgsmrange(2)-Zgsmrange(1)) )

ip=(ipx-1)*YLEN*ZLEN+(ipy-1)*ZLEN+ipz-1;

xmin=double(0+(ipx-1)*Lx/XLEN)
xmax=double(xmin+Lx/XLEN)
ymin=double(0+(ipy-1)*Ly/YLEN)
ymax=double(ymin+double(Ly/YLEN));
zmin=double(0+(ipz-1)*Lz/ZLEN)
zmax=double(zmin+double(Lz/ZLEN))

nx=Nxc/XLEN+1;
ny=Nyc/YLEN+1;
nz=Nzc/ZLEN+1;

dx=double(xmax-xmin)/double(nx-1);
dy=double(ymax-ymin)/double(ny-1);
dz=double(zmax-zmin)/double(nz-1);

[xg,yg,zg]=ndgrid(xmin:dx:xmax,ymin:dy:ymax,zmin:dz:zmax);



for is=1:1

if(mod(is,2)==0)
%vmax = .03/3;
vmax=.07*3;
vmin = -vmax;
else
%vmax = .2/3;
vmax=.4
vmin = -vmax;
end

partfile='restart';


info=hdf5info([results_dir partfile num2str(ip) '.hdf']);
[nt dummy]=size(info.GroupHierarchy.Groups(3).Groups(is).Groups(5).Datasets(:));


it=nt;
%ex=squeeze(hdf5read(infof.GroupHierarchy.Groups(1).Groups(4).Datasets(it)));


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

david=1;
if(david)
%David's coordinates
bp2D=sqrt(bxp.^2+byp.^2);
perp2x=bzp.*bxp./(bp.*bp2D);
perp2y=bzp.*byp./(bp.*bp2D);
perp2z=-bp2D./bp;
else
% max p pane coordinates
bp2D=sqrt(bxp.^2+byp.^2);
perp1x=bzp./sqrt(bzp.^2+bxp.^2);
perp1z=-bxp./sqrt(bzp.^2+bxp.^2);
perp2x=-bxp.*byp;
perp2y=(bxp.^2+bzp.^2);
perp2z=-bzp.*byp;
denom=sqrt(perp2x.^2+perp2y.^2+perp2z.^2);
perp2x=perp2x./denom;
perp2y=perp2y./denom;
perp2z=perp2z./denom;
end

u=squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is).Groups(2).Datasets(it)));
v=squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is).Groups(3).Datasets(it)));
w=squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is).Groups(4).Datasets(it)));

upar=(u.*bxp+v.*byp+w.*bzp)./(bp+1e-10);

if(david)
%David's coordinates
uperp1=(byp.*u-bxp.*v)./bp2D;
uperp2=perp2x.*u+perp2y.*v+perp2z.*w;
else
uperp1=perp1x.*u+perp1z.*w;
uperp2=perp2x.*u+perp2y.*v+perp2z.*w;
end
%
%	Subdivide the processors
%
nsubx=8;
nsuby=4;
nsubz=1;
dxsub=(xmax-xmin)/nsubx;
dysub=(ymax-ymin)/nsuby;
dzsub=(zmax-zmin)/nsubz;
ixsub=floor((x-xmin)/dxsub);
iysub=floor((y-ymin)/dysub);
izsub=floor((z-zmin)/dzsub);

% choose the desired subdomain in y ans z, from 0 to nsub-1
suby=0;subz=0;
for subx=0:nsubx-1
ii=(ixsub==subx)&(iysub==suby)&(izsub==subz);

bxscan((ipx-1)*nsubx+subx+1)=mean(bxp(ii));
byscan((ipx-1)*nsubx+subx+1)=mean(byp(ii));
bzscan((ipx-1)*nsubx+subx+1)=mean(bzp(ii));

xscan((ipx-1)*nsubx+subx+1)=(mean(x(ii)));
yscan((ipx-1)*nsubx+subx+1)=(mean(y(ii)));
zscan((ipx-1)*nsubx+subx+1)=(mean(z(ii)));


%xscan((ipx-1)*nsub+subx+1)=gsmx(mean(x(ii)));
%yscan((ipx-1)*nsub+subx+1)=gsmx(mean(y(ii)));
%zscan((ipx-1)*nsub+subx+1)=gsmx(mean(z(ii)));

ndiv=101;
vdf_sp=spaziofasi3D(upar(ii),uperp1(ii),uperp2(ii),q(ii),vmin,vmax,ndiv);
vdf_sp=vdf_sp./sum(vdf_sp(:));

% converting to physical units
tail

Qp=abs(sum(q(ii)));
rho=Qp./(xmax-xmin)/(ymax-ymin)/(zmax-zmin);

vdf_sp=vdf_sp*rho*np/c^3;

%filename=['THOR_' 'species_' num2str(is) '_' num2str(ipx*nsub+subx) '.txt'];
%saveTHOR(vdf_sp, xscan((ipx-1)*nsub+subx+1)*dp, filename, ...
%         vmin*c,vmin*c,vmin*c,(vmax-vmin)*c,(vmax-vmin)*c,(vmax-vmin)*c, ...
%         bxscan((ipx-1)*nsub+subx+1)*Bnorm,byscan((ipx-1)*nsub+subx+1)*Bnorm,bzscan((ipx-1)*nsub+subx+1)*Bnorm);


vdf_sp=smooth3D(vdf_sp,3);
vdf_sp=vdf_sp./sum(vdf_sp(:));
dv=(vmax-vmin)/ndiv;
%filename=['vdf_' 'species_' num2str(is) '_' num2str(ipx*nsub+subx) '.vtk'];
%savevtk_bin(vdf_sp,filename,'vdf',dv,dv,dv,vmin,vmin,vmin);





global labelT
labelT='x/R_E=';
labelT='x/d_i=';
symmetric_color=0
color_choice=4

% for 101x101 plots
if(is==1)
crange=[-8 -2];
else
crange=[-8 -2];
end 
crange=[-6 -2];

immagine_dir([vmin vmax],[vmin vmax],log10(1e-20+squeeze(sum(vdf_sp,3))), ...
             ['vdfParPerp_' 'species_' num2str(is) '_' num2str(ipx*nsubx+subx)], ...
             crange,0,num2str(xscan((ipx-1)*nsubx+subx+1)),max(size(q(ii))),1,'v_{||}/c','v_{\perp 1}/c','vdf');

%immagine_dir([vmin vmax],[vmin vmax],(1e-10+squeeze(sum(vdf_sp,2))+squeeze(sum(vdf_sp,3)))/2, ...
%             ['vdfXZ_' 'species_' num2str(is) '_' num2str(ipx*nsub+subx)], ...
%             [0 0],0,num2str(xscan((ipx-1)*nsub+subx+1)),max(size(q(ii))),1,'v_{||}/c','v_{\perp,avg}/c','vdf');

immagine_dir([vmin vmax],[vmin vmax],log10(1e-20+squeeze(sum(vdf_sp,1))), ...
             ['vdfPerp1Perp2_' 'species_' num2str(is) '_' num2str(ipx*nsubx+subx)], ...
             crange,0,num2str(xscan((ipx-1)*nsubx+subx+1)),max(size(q(ii))),2,'v_{\perp 1}/c','v_{\perp 2}/c','vdf');
end


%plot(xscan,bxscan,xscan,byscan,xscan,bzscan)
%legend('Bx','By','Bz') 
%print -dpng scanni.png
end

end
