function [Pxx,Pxy,Pxz,Pyy,Pyz,Pzz,Ppar,Pper1,Pper2,Peps,Nx,Ny,Nz] = read_binVTK_pressure(dir,name,cycle)
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
s = fgetl(fid); % LOOKUP_TABLE default

Pxx=fread(fid,Nx*Ny*Nz,'float','b'); %'b' means big endian
Pxx=reshape(Pxx,Nx,Ny,Nz);

s = fgetl(fid); % SCALARS/VECTORS name data_type (ex: SCALARS imagedata unsigned_char)
s = fgetl(fid); % LOOKUP_TABLE default

Pxy=fread(fid,Nx*Ny*Nz,'float','b'); %'b' means big endian
Pxy=reshape(Pxy,Nx,Ny,Nz);

s = fgetl(fid); % SCALARS/VECTORS name data_type (ex: SCALARS imagedata unsigned_char)
s = fgetl(fid); % LOOKUP_TABLE default

Pxz=fread(fid,Nx*Ny*Nz,'float','b'); %'b' means big endian
Pxz=reshape(Pxz,Nx,Ny,Nz);

s = fgetl(fid); % SCALARS/VECTORS name data_type (ex: SCALARS imagedata unsigned_char)
s = fgetl(fid); % LOOKUP_TABLE default

Pyy=fread(fid,Nx*Ny*Nz,'float','b'); %'b' means big endian
Pyy=reshape(Pyy,Nx,Ny,Nz);

s = fgetl(fid); % SCALARS/VECTORS name data_type (ex: SCALARS imagedata unsigned_char)
s = fgetl(fid); % LOOKUP_TABLE default

Pyz=fread(fid,Nx*Ny*Nz,'float','b'); %'b' means big endian
Pyz=reshape(Pyz,Nx,Ny,Nz);

s = fgetl(fid); % SCALARS/VECTORS name data_type (ex: SCALARS imagedata unsigned_char)
s = fgetl(fid); % LOOKUP_TABLE default

Pzz=fread(fid,Nx*Ny*Nz,'float','b'); %'b' means big endian
Pzz=reshape(Pzz,Nx,Ny,Nz);

s = fgetl(fid); % SCALARS/VECTORS name data_type (ex: SCALARS imagedata unsigned_char)
s = fgetl(fid); % LOOKUP_TABLE default

Ppar=fread(fid,Nx*Ny*Nz,'float','b'); %'b' means big endian
Ppar=reshape(Ppar,Nx,Ny,Nz);

s = fgetl(fid); % SCALARS/VECTORS name data_type (ex: SCALARS imagedata unsigned_char)
s = fgetl(fid); % LOOKUP_TABLE default

Pper1=fread(fid,Nx*Ny*Nz,'float','b'); %'b' means big endian
Pper1=reshape(Pper1,Nx,Ny,Nz);

s = fgetl(fid); % SCALARS/VECTORS name data_type (ex: SCALARS imagedata unsigned_char)
s = fgetl(fid); % LOOKUP_TABLE default

Pper2=fread(fid,Nx*Ny*Nz,'float','b'); %'b' means big endian
Pper2=reshape(Pper2,Nx,Ny,Nz);

s = fgetl(fid); % SCALARS/VECTORS name data_type (ex: SCALARS imagedata unsigned_char)
s = fgetl(fid); % LOOKUP_TABLE default

Peps=fread(fid,Nx*Ny*Nz,'float','b'); %'b' means big endian
Peps=reshape(Peps,Nx,Ny,Nz);
fclose(fid);
end

