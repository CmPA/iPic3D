function [Vx, Vy, Vz,Nx,Ny,Nz] = read_binVTK_vector(dir,name,cycle)
%Reads the binary VTKs

ncycle=num2str(cycle);

filename=[dir name '_xyz_cycle' ncycle '.vtk']
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


V=fread(fid,3*Nx*Ny*Nz,'float','b'); %'b' means big endian

Vx=reshape(V(1:3:end),Nx,Ny,Nz);
Vy=reshape(V(2:3:end),Nx,Ny,Nz);
Vz=reshape(V(3:3:end),Nx,Ny,Nz);
fclose(fid);
end

