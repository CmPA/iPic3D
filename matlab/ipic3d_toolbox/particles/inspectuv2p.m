addpath '/home/gianni/matlab2/matlab-parsek'
addpath '/home/gianni/matlab2'
!rm cazzo
close all

!h5dump  -o cazzo -b -y -d /particles/species_2/q/cycle_80000 /shared/gianni/run100/part975.hdf 
!h5dump  -o cazzo2 -b -y -d /particles/species_2/q/cycle_80000 /shared/gianni/run100/part976.hdf 

fid=fopen('cazzo','r');
q1=fread(fid,inf,'float64');
fclose(fid)

fid=fopen('cazzo2','r');
q2=fread(fid,inf,'float64');
fclose(fid)

!h5dump  -o cazzo -b -y -d /particles/species_0/q/cycle_80000 /shared/gianni/run100/part975.hdf 
!h5dump  -o cazzo2 -b -y -d /particles/species_0/q/cycle_80000 /shared/gianni/run100/part976.hdf 

fid=fopen('cazzo','r');
q3=fread(fid,inf,'float64');
fclose(fid)

fid=fopen('cazzo2','r');
q4=fread(fid,inf,'float64');
fclose(fid)

q=[q1;q2;q3;q4];

!h5dump  -o cazzo -b -y -d /particles/species_2/x/cycle_80000 /shared/gianni/run100/part975.hdf 
!h5dump  -o cazzo2 -b -y -d /particles/species_2/x/cycle_80000 /shared/gianni/run100/part976.hdf 

fid=fopen('cazzo','r');
x1=fread(fid,inf,'float64');
fclose(fid)

fid=fopen('cazzo2','r');
x2=fread(fid,inf,'float64');
fclose(fid)

!h5dump  -o cazzo -b -y -d /particles/species_0/x/cycle_80000 /shared/gianni/run100/part975.hdf 
!h5dump  -o cazzo2 -b -y -d /particles/species_0/x/cycle_80000 /shared/gianni/run100/part976.hdf 

fid=fopen('cazzo','r');
x3=fread(fid,inf,'float64');
fclose(fid)

fid=fopen('cazzo2','r');
x4=fread(fid,inf,'float64');
fclose(fid)

x=[x1;x2;x3;x4];

!h5dump  -o cazzo -b -y -d /particles/species_2/y/cycle_80000 /shared/gianni/run100/part975.hdf 
!h5dump  -o cazzo2 -b -y -d /particles/species_2/y/cycle_80000 /shared/gianni/run100/part976.hdf 

fid=fopen('cazzo','r');
x1=fread(fid,inf,'float64');
fclose(fid)

fid=fopen('cazzo2','r');
x2=fread(fid,inf,'float64');
fclose(fid)

!h5dump  -o cazzo -b -y -d /particles/species_0/y/cycle_80000 /shared/gianni/run100/part975.hdf 
!h5dump  -o cazzo2 -b -y -d /particles/species_0/y/cycle_80000 /shared/gianni/run100/part976.hdf 

fid=fopen('cazzo','r');
x3=fread(fid,inf,'float64');
fclose(fid)

fid=fopen('cazzo2','r');
x4=fread(fid,inf,'float64');
fclose(fid)

y=[x1;x2;x3;x4];

!h5dump  -o cazzo -b -y -d /particles/species_2/w/cycle_80000 /shared/gianni/run100/part975.hdf 
!h5dump  -o cazzo2 -b -y -d /particles/species_2/w/cycle_80000 /shared/gianni/run100/part976.hdf

fid=fopen('cazzo','r');
u1=fread(fid,inf,'float64');
fclose(fid)

fid=fopen('cazzo2','r');
u2=fread(fid,inf,'float64');
fclose(fid)

!h5dump  -o cazzo -b -y -d /particles/species_0/w/cycle_80000 /shared/gianni/run100/part975.hdf 
!h5dump  -o cazzo2 -b -y -d /particles/species_0/u/cycle_80000 /shared/gianni/run100/part976.hdf

fid=fopen('cazzo','r');
u3=fread(fid,inf,'float64');
fclose(fid)

fid=fopen('cazzo2','r');
u4=fread(fid,inf,'float64');
fclose(fid)


u=[u1;u2;u3;u4];

ii=y<10.1 & y>9.9;
x=x(ii);y=y(ii);u=u(ii);q=q(ii);

figure(1)
plot(x,u,'.','MarkerSize',1)
axis equal
%figure(2)
%hist(u,100)

figure(2)
[totnum,nbinu,xrange,urange]=spaziofasi(x',u',q',0);
surf(xrange(1:end-1)',urange(1:end-1)',nbinu(1:end-1,1:end-1)','edgecolor','none','facecolor','blue')
lighting phong
shading interp
camlight(0,90) % luce dall'alto
view(2) %visione dall'alto
axis tight
