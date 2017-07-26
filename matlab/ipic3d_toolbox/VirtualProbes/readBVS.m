fid=fopen('VirtualSatelliteTraces0.txt')
ennes=fscanf(fid,'%f',4) 
nsat=ennes(1) 
nsatx=ennes(2); 
nsaty=ennes(3);
nsatz=ennes(4);
for i=1:nsat
x=fscanf(fid,'%f',3); 
xp(i)=x(1);
yp(i)=x(2);
zp(i)=x(3);
end
fscanf(fid,'%f');
a=fread(fid,'float');
fclose(fid)
