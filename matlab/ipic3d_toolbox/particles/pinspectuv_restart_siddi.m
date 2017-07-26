addpath '/home/gianni/matlab2/matlab-parsek'
addpath '/home/gianni/matlab2'

% This works i the run has saved the particles
%results_dir='/shared/gianni/lhdi6redux/';
%results_dir='/shared02/gianni/tred60/data/';
%results_dir='/shared/gianni/laila11b/';
results_dir = '/home/gianni/';

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

ipx=XLEN/2;
%for  ipx=1:XLEN

% Per siddi
for ipx=1:1
%ipy=YLEN/2+2;
ipy=YLEN/2+2;
ipz=ZLEN/2;

% Per Siddi
ipy=1;
ipz=1;
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

[zg,yg,xg]=ndgrid(zmin:dz:zmax,ymin:dy:ymax,xmin:dx:xmax);


close all

vmax_graph=[.25 .06 .25 .06];
for is=1:2

vmax = vmax_graph(is);
partfile='part';
partfile='restart';
if exist([results_dir partfile num2str(ip) '.hdf']) ==2 

info=hdf5info([results_dir partfile num2str(ip) '.hdf']);
[nt dummy]=size(info.GroupHierarchy.Groups(3).Groups(is).Groups(5).Datasets(:));


infof=hdf5info([results_dir 'proc' num2str(ip) '.hdf']);

[ntf dummy]=size(infof.GroupHierarchy.Groups(1).Groups(4).Datasets(:));

for it=nt:nt
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
bxp=interpn(xg,yg,zg,Bx,x,y,z);
byp=interpn(xg,yg,zg,By,x,y,z);
bzp=interpn(xg,yg,zg,Bz,x,y,z);
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


if (is == 1) 
histpar=[histpar; hist(epar,bins)];
histperp=[histperp; hist(eperp,bins)]; 
xhist=[xhist;mean(x)];
bhist=[bhist;mean(bp)];
bxhist=[bxhist;mean(bxp)];
byhist=[byhist;mean(byp)];
bzhist=[bzhist;mean(bzp)];
end

max(y)
min(y)

%ii=y<7;
%q=q(ii);
%x=x(ii);
%y=y(ii);
%u=u(ii);
% q2=squeeze(hdf5read(info2.GroupHierarchy.Groups(1).Groups(1).Groups(2).Datasets(it)));
% x2=squeeze(hdf5read(info2.GroupHierarchy.Groups(1).Groups(1).Groups(6).Datasets(it)));
% u2=squeeze(hdf5read(info2.GroupHierarchy.Groups(1).Groups(1).Groups(3).Datasets(it)));
% q=[q;q2];
% x=[x;x2];
% u=[u;u2];
%figure(1)
%plot(x,y,'.','MarkerSize',1)
%axis square
%set(gcf,'Renderer','zbuffer');
%print('-dpng',['filename '.png'])
%print -dpng figure1


figure(1)
[totnum,nbinu,xrange,urange]=spaziofasi(upar',uperp1',q',0);
surf(xrange(1:end-1)',urange',log(nbinu(1:end-1,:))','edgecolor','none','facecolor','blue')
lighting phong
shading interp
camlight(0,90) % luce dall'alto
view(2) %visione dall'alto
xlim([-vmax vmax])
ylim([-vmax vmax])	 
	 axis square 
title(['x=' num2str((xmin+xmax)/2) '  y=' num2str((ymin+ymax)/2) '  z=' num2str((zmin+zmax)/2) ])% '  xp=' num2str(mean(x))  '   yp=' num2str(mean(y))]) 
xlabel('upar')
	 ylabel('uperp1')
	 colorbar
set(gcf,'Renderer','zbuffer');
print('-dpng',['SP' num2str(is) 'phasespaceUV' num2str(ipx) '.png'])


figure(1)
[totnum,nbinu,xrange,urange]=spaziofasi(upar',uperp2',q',0);
surf(xrange(1:end-1)',urange',log(nbinu(1:end-1,:))','edgecolor','none','facecolor','blue')
lighting phong
shading interp
camlight(0,90) % luce dall'alto
view(2) %visione dall'alto
	 xlim([-vmax vmax])
	 ylim([-vmax vmax])	 
	 axis square 
title(['x=' num2str((xmin+xmax)/2) '  y=' num2str((ymin+ymax)/2) '  z=' num2str((zmin+zmax)/2) ])% '  xp=' num2str(mean(x))  '   yp=' num2str(mean(y))]) 
	 xlabel('upar')
	 ylabel('uperp2')
	 colorbar
set(gcf,'Renderer','zbuffer');
print('-dpng',['SP' num2str(is) 'phasespaceUW' num2str(ipx) '.png'])

	 figure(1)
	 [totnum,nbinu,xrange,urange]=spaziofasi(uperp1',uperp2',q',0);
				surf(xrange(1:end-1)',urange',log(nbinu(1:end-1,:))','edgecolor','none','facecolor','blue')
					 lighting phong
					 shading interp
					 camlight(0,90) % luce dall'alto
					 view(2) %visione dall'alto
					 xlim([-vmax vmax])
					 ylim([-vmax vmax])	 
					 axis square 
					 title(['x=' num2str((xmin+xmax)/2) '  y=' num2str((ymin+ymax)/2) '  z=' num2str((zmin+zmax)/2) ])% '  xp=' num2str(mean(x))  '   yp=' num2str(mean(y))]) 
					 xlabel('uperp1')
					 ylabel('uperp2')
					 colorbar
					 set(gcf,'Renderer','zbuffer');
					 print('-dpng',['SP' num2str(is) 'phasespaceUW' num2str(ipx) '.png'])
					 

% This below works only in 2D
%[nx ny]=size(ex);
%dx=double(xmax-xmin)/double(nx-1);
%dy=double(ymax-ymin)/double(ny-1);
%[xg,yg]=meshgrid(xmin:dx:xmax,ymin:dy:ymax);

%figure(3)
%pcolor(xg,yg,ex)
%shading interp
%set(gcf,'Renderer','zbuffer');
%print -dpng figure3
%pause
end

end

end

end

figure
subplot(2,1,1)
for i=1:1
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
for i=1:1
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


ii=[1 5 6 8];

figure
semilogy(bins,histpar(ii,:))
legend(num2str(xhist(ii)))
title('Parallel')
xlabel('Energy')
ylabel('pdf(E)') 
set(gcf,'Renderer','zbuffer');
print('-dpng','tmp.png')


figure
alfa= flipud(cumsum(flipud(histpar')))';
for i=1:1
alfa(i,:)=alfa(i,:)/alfa(i,1);
end

loglog(bins,alfa(ii,:))
xlim([1e-4 1e-1])
legend(num2str(xhist(ii))) 
title('Parallel')
xlabel('Energy')
ylabel('cdf(E)') 
set(gcf,'Renderer','zbuffer');
print('-dpng','tmp2.png')

figure
semilogy(bins,histperp(ii,:))
legend(num2str(xhist(ii))) 
set(gcf,'Renderer','zbuffer');
print('-dpng','tmp3.png')


figure
alfa= flipud(cumsum(flipud(histperp')))';
for i=1:1
alfa(i,:)=alfa(i,:)/alfa(i,1);
end
ii=[1 5 6 7];
loglog(bins,alfa(ii,:),bins.*byhist(5)./byhist(7),alfa(7,:),'--',bins,alfa(7,:),'x')
xlim([1e-4 1e-1])
legend(num2str(xhist(ii))) 
set(gcf,'Renderer','zbuffer');
print('-dpng','tmp4.png')
