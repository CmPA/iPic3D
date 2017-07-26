function [V,Nx,Ny,Nz,dx,dy,dz] = read_binVTK_scalar(dir,name,cycle)
%Reads the binary VTKs

ncycle=num2str(cycle);

%filename=[dir name '_xyz_cycle' ncycle '.vtk']
filename=[dir name  ncycle '.vtk']
fid = fopen(filename,'r');
fgetl(fid); % # vtk DataFile Version x.x
fgetl(fid); % comments
fgetl(fid); % ASCII
fgetl(fid); % DATASET STRUCTURED_POINTS

s = fgetl(fid); % DIMENSIONS NX NY NZ
sz = sscanf(s, '%*s%d%d%d');
Nx=sz(1);
Ny=sz(2);
Nz=sz(3);

s=fgetl(fid); % ORIGIN OX OY OZ
oo = sscanf(s, '%*s%d%d%d');
s=fgetl(fid); % SPACING SX SY SZ
dd= sscanf(s, '%*s%f%f%f');
dx=dd(1);
dy=dd(2);
dz=dd(3);
s=fgetl(fid); % POINT_DATA NXNYNZ
npoints = sscanf(s, '%*s%d%d%d');

s = fgetl(fid); % SCALARS/VECTORS name data_type (ex: SCALARS imagedata unsigned_char)
s = fgetl(fid); % LOOKUP_TABLE default

V=fread(fid,Nx*Ny*Nz,'float','b'); %'b' means big endian
V=reshape(V,Nx,Ny,Nz);
fclose(fid);
end

