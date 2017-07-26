function [V1,V2,Nx,Ny,Nz] = read_binVTK_multiscalar(dir,name,cycle)
%Reads the binary VTKs

ncycle=num2str(cycle);

filename=[dir name '_xyz_cycle' ncycle '.vtk'];
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
s=fgetl(fid); % POINT_DATA NXNYNZ
npoints = sscanf(s, '%*s%d%d%d');

s = fgetl(fid); % SCALARS/VECTORS name data_type (ex: SCALARS imagedata unsigned_char)
svstr = sscanf(s, '%s', 1);
dtstr = sscanf(s, '%*s%*s%s');
s = fgetl(fid); % LOOKUP_TABLE default

V1=fread(fid,Nx*Ny*Nz,'float','b'); %'b' means big endian
V1=reshape(V1,Nx,Ny,Nz);

s = fgetl(fid); % SCALARS/VECTORS name data_type
s = fgetl(fid); % LOOKUP_TABLE default

V2=fread(fid,Nx*Ny*Nz,'float','b'); %'b' means big endian
V2=reshape(V2,Nx,Ny,Nz);

fclose(fid);
end

