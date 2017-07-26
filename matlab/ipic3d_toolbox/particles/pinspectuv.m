addpath '/home/gianni/matlab2/matlab-parsek'
addpath '/home/gianni/matlab2'

% This works i the run has saved the particles
%results_dir='/shared/gianni/lhdi6redux/';
results_dir='/shared02/gianni/tred60/data/';

filename=[results_dir 'settings.hdf'];
XLEN=double(hdf5read(filename,'/topology/XLEN'))
YLEN=double(hdf5read(filename,'/topology/YLEN'))
Lx=double(hdf5read(filename,'/collective/Lx'))
Ly=double(hdf5read(filename,'/collective/Ly'))
B0x=double(hdf5read(filename,'/collective/Bx0'))

% use for pqrsek 2D 
%B0x=double(hdf5read(filename,'/collective/B0x'))

Dt=double(hdf5read(filename,'/collective/Dt'))
ipx=XLEN/2;
ipy=YLEN/2+2;
ip=(ipx-1)*YLEN+(ipy-1)
ipx2=ipx;
ipy2=ipy-1;
ip2=(ipx2-1)*YLEN+ipy2

xmin=double(0+(ipx-1)*Lx/XLEN)
xmax=double(xmin+Lx/XLEN)
ymin=double(0+(ipy-1)*Ly/YLEN)
ymax=double(ymin+double(Ly/YLEN))

close all

is=3
partfile='part';
partfile='restart';
info=hdf5info([results_dir partfile num2str(ip) '.hdf']);
[nt dummy]=size(info.GroupHierarchy.Groups(1).Groups(is).Groups(5).Datasets(:));
info2=hdf5info([results_dir partfile num2str(ip2) '.hdf']);


infof=hdf5info([results_dir 'proc' num2str(ip) '.hdf']);

[ntf dummy]=size(infof.GroupHierarchy.Groups(2).Groups(4).Datasets(:));

for it=nt:nt
ex=squeeze(hdf5read(infof.GroupHierarchy.Groups(2).Groups(4).Datasets(it)));   

q=squeeze(hdf5read(info.GroupHierarchy.Groups(1).Groups(is).Groups(2).Datasets(it)));
x=squeeze(hdf5read(info.GroupHierarchy.Groups(1).Groups(is).Groups(6).Datasets(it)));
y=squeeze(hdf5read(info.GroupHierarchy.Groups(1).Groups(is).Groups(7).Datasets(it)));
u=squeeze(hdf5read(info.GroupHierarchy.Groups(1).Groups(is).Groups(4).Datasets(it)));

max(y)
min(y)

ii=y<7;
q=q(ii);
x=x(ii);
y=y(ii);
u=u(ii);
% q2=squeeze(hdf5read(info2.GroupHierarchy.Groups(1).Groups(1).Groups(2).Datasets(it)));
% x2=squeeze(hdf5read(info2.GroupHierarchy.Groups(1).Groups(1).Groups(6).Datasets(it)));
% u2=squeeze(hdf5read(info2.GroupHierarchy.Groups(1).Groups(1).Groups(3).Datasets(it)));
% q=[q;q2];
% x=[x;x2];
% u=[u;u2];
figure(1)
plot(x,y,'.','MarkerSize',1)
axis equal
set(gcf,'Renderer','zbuffer');
%print('-dpng',[filename '.png'])
print -dpng figure1


figure(2)
[totnum,nbinu,xrange,urange]=spaziofasi(y',u',q',0);
surf(xrange',urange',nbinu','edgecolor','none','facecolor','blue')
lighting phong
shading interp
camlight(0,90) % luce dall'alto
view(2) %visione dall'alto
axis tight
set(gcf,'Renderer','zbuffer');
print -dpng figure2


[nx ny]=size(ex);
dx=double(xmax-xmin)/double(nx-1);
dy=double(ymax-ymin)/double(ny-1);
[xg,yg]=meshgrid(xmin:dx:xmax,ymin:dy:ymax);

figure(3)
pcolor(xg,yg,ex)
shading interp
set(gcf,'Renderer','zbuffer');
print -dpng figure3
%pause
end
