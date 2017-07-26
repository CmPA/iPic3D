!rm cazzo
close all
!h5dump  -o cazzo -b -y -d /particles/species_0/x/cycle_60000 /shared/gianni/run100/part1009.hdf 

fid=fopen('cazzo','r');
x=fread(fid,inf,'float64');
fclose(fid)

!h5dump  -o cazzo -b -y -d /particles/species_0/y/cycle_60000 /shared/gianni/run100/part1009.hdf 

fid=fopen('cazzo','r');
u=fread(fid,inf,'float64');
fclose(fid)

figure(1)
plot(x,u,'.','MarkerSize',1)
axis equal
%figure(2)
%hist(u,100)

figure(2)
[totnum,nbinu,xrange,urange]=spaziofasi(x',u',ones(size(x')),0);
surf(xrange',urange',nbinu','edgecolor','none','facecolor','blue')
lighting phong
shading interp
camlight(0,90) % luce dall'alto
view(2) %visione dall'alto
axis tight
