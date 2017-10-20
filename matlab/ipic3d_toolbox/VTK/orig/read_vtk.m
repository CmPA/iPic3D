function [V,Vx,Vy,Vz,dx,dy,dz]=read_vtk(filename)
V=[];
Vx=[];
Vy=[];
Vz=[];
fid = fopen(filename,'r');
fgetl(fid); % # vtk DataFile Version x.x
fgetl(fid); % comments
fgetl(fid); % ASCII
fgetl(fid); % DATASET STRUCTURED_POINTS

s = fgetl(fid); % DIMENSIONS NX NY NZ
sz = sscanf(s, '%*s%d%d%d')

s=fgetl(fid); % ORIGIN OX OY OZ
oo = sscanf(s, '%*s%d%d%d')
s=fgetl(fid); % SPACING SX SY SZ
dd= sscanf(s, '%*s%f%f%f')
s=fgetl(fid); % POINT_DATA NXNYNZ
npoints = sscanf(s, '%*s%d%d%d')

s = fgetl(fid); % SCALARS/VECTORS name data_type (ex: SCALARS imagedata unsigned_char)
svstr = sscanf(s, '%s', 1)
dtstr = sscanf(s, '%*s%*s%s')

if( strcmp(svstr,'SCALARS') > 0 )
fgetl(fid); % the lookup table
V=[];
for np=1:npoints
s=fgetl(fid); 
V = [V;sscanf(s, '%g').'];
end
V=reshape(V,sz(1),sz(2));

elseif( strcmp(svstr,'VECTORS') > 0 )
% read data
V=[];
for np=1:npoints
s=fgetl(fid); 
V = [V;sscanf(s, '%g%g%g').'];
end
     Vx=reshape(V(:,1),sz(1),sz(2));
     Vy=reshape(V(:,2),sz(1),sz(2));
     Vz=reshape(V(:,3),sz(1),sz(2));
     V=sqrt(Vx.^2+Vy.^2+Vz.^2);
end

V=V';

    Vx=Vx';
    Vy=Vy';
    Vz=Vz';
   
dx=dd(1);
dy=dd(2);
dz=dd(3);
fclose(fid);



