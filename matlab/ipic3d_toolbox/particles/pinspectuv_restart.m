addpath '/home/gianni/matlab2/matlab-parsek'
addpath '/home/gianni/matlab2'

% This works i the run has saved the particles
%results_dir='/shared/gianni/lhdi6redux/';
results_dir='/shared02/gianni/tred60/data/';

filename=[results_dir 'settings.hdf'];
XLEN=double(hdf5read(filename,'/topology/XLEN'))
YLEN=double(hdf5read(filename,'/topology/YLEN'))
ZLEN=double(hdf5read(filename,'/topology/ZLEN'))
Lx=double(hdf5read(filename,'/collective/Lx'))
Ly=double(hdf5read(filename,'/collective/Ly'))
Lz=double(hdf5read(filename,'/collective/Lz'))
B0x=double(hdf5read(filename,'/collective/Bx0'))

% use for pqrsek 2D 
%B0x=double(hdf5read(filename,'/collective/B0x'))

Dt=double(hdf5read(filename,'/collective/Dt'))
ipx=XLEN/2;
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

close all

for is=1:4

partfile='part';
partfile='restart';
if exist([results_dir partfile num2str(ip) '.hdf']) ==2 

info=hdf5info([results_dir partfile num2str(ip) '.hdf']);
[nt dummy]=size(info.GroupHierarchy.Groups(3).Groups(is).Groups(5).Datasets(:));


infof=hdf5info([results_dir 'proc' num2str(ip) '.hdf']);

[ntf dummy]=size(infof.GroupHierarchy.Groups(1).Groups(4).Datasets(:));

for it=nt:nt
ex=squeeze(hdf5read(infof.GroupHierarchy.Groups(1).Groups(4).Datasets(it)));   

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
u=squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is).Groups(2).Datasets(it)));
v=squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is).Groups(3).Datasets(it)));
w=squeeze(hdf5read(info.GroupHierarchy.Groups(3).Groups(is).Groups(4).Datasets(it)));

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
%axis equal
%set(gcf,'Renderer','zbuffer');
%print('-dpng',['filename '.png'])
%print -dpng figure1


figure(1)
[totnum,nbinu,xrange,urange]=spaziofasi(u',v',q',0);
surf(xrange(1:end-1)',urange',log(nbinu(1:end-1,:))','edgecolor','none','facecolor','blue')
lighting phong
shading interp
camlight(0,90) % luce dall'alto
view(2) %visione dall'alto
axis tight
title(['x=' num2str((xmin+xmax)/2) '  y=' num2str((ymin+ymax)/2) '  z=' num2str((zmin+zmax)/2) ])% '  xp=' num2str(mean(x))  '   yp=' num2str(mean(y))]) 
colorbar
set(gcf,'Renderer','zbuffer');
print('-dpng',['SP' num2str(is) 'phasespaceUV' num2str(ipx) '.png'])


figure(1)
[totnum,nbinu,xrange,urange]=spaziofasi(u',w',q',0);
surf(xrange(1:end-1)',urange',log(nbinu(1:end-1,:))','edgecolor','none','facecolor','blue')
lighting phong
shading interp
camlight(0,90) % luce dall'alto
view(2) %visione dall'alto
axis tight
title(['x=' num2str((xmin+xmax)/2) '  y=' num2str((ymin+ymax)/2) '  z=' num2str((zmin+zmax)/2) ])% '  xp=' num2str(mean(x))  '   yp=' num2str(mean(y))]) 
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
