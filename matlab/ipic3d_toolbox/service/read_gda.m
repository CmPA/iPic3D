function [V]=read_gda(dir,name,cycle,Nx,Ny,Nz)

ncycle=num2str(cycle);
file=[dir name '_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
V=fread(fid,'real*8');
fclose(fid);
V=reshape(V,Nx,Ny,Nz);


