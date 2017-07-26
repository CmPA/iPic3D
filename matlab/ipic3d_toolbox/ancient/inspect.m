!h5dump  -o cazzo -b -y -d /particles/species_0/x/cycle_0 proc0.hdf 

fid=fopen('cazzo','r');
x=fread(fid,inf,'float64');
fclose(fid)

!h5dump  -o cazzo -b -y -d /particles/species_0/w/cycle_0 proc0.hdf 

fid=fopen('cazzo','r');
u=fread(fid,inf,'float64');
fclose(fid)

figure(1)
plot(x,u,'.','MarkerSize',1)
figure(2)
hist(u,100)
