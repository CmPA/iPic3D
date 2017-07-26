
for cycle=Ncyc_ini:1000:Ncyc_max

ncycle=num2str(cycle)
ncycle1=num2str(cycle,'%06d')

file=[dir 'B_x_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Bx=fread(fid,'real*8');
fclose(fid);
Bx=reshape(Bx,Nx,Ny,Nz);

file=[dir 'B_y_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
By=fread(fid,'real*8');
fclose(fid);
By=reshape(By,Nx,Ny,Nz);


file=[dir 'B_z_cycle' ncycle '.gda'];
fid= fopen(file,'rb');
Bz=fread(fid,'real*8');
fclose(fid);
Bz=reshape(Bz,Nx,Ny,Nz);

Nsm=3
Bxsm=smooth3D(Bx,Nsm);
Bysm=smooth3D(By,Nsm);
Bzsm=smooth3D(Bz,Nsm);

ox=max(Xgsmrange);
oy=min(Zgsmrange);
oz=min(Ygsmrange);
ddx=dx/Lx*(Xgsmrange(2)-Xgsmrange(1));
ddy=dy/Ly*(Zgsmrange(2)-Zgsmrange(1));
ddz=dz/Lz*(Ygsmrange(2)-Ygsmrange(1));

savevtkvector(Bx,  By, Bz, ['B' ncycle1 '.vtk'],'B',ddx,ddy,ddz,ox,oy,oz)

end