
addpath(genpath('Users/gianni/Dropbox/Science/codes/matlab4iPic3D/ipic3d_toolbox'));

clear all
close all


yesn=0
if(yesn)
file='/Volumes/Granata/paraview/tred54/rho_xyz_cycle15000.vtk'

[ne,ni]=read_vtk_multiscalar_3d(file,2);

[Nx,Ny,Nz]=size(ne);

ne=circshift(ne,Nx/2,1);
ni=circshift(ni,Nx/2,1);

dx= 0.078125;
dy= 0.078125;
dz=0.078125;

savevtk_bin(ni, 'Ni15000.vtk','ni',dx,dy,dz,0,0,0)
savevtk_bin(ne, 'Ne15000.vtk','ne',dx,dy,dz,0,0,0)
end

yesB=0
if(yesB)
file='/Volumes/Granata/paraview/tred54/B_xyz_cycle15000.vtk'

[V,Bx,By,Bz]=read_vtk_3d(file,0);

[Nx,Ny,Nz]=size(Bx);

Vx=circshift(Bx,Nx/2,1);
Vy=circshift(By,Nx/2,1);
Vz=circshift(Bz,Nx/2,1);

dx= 0.078125;
dy= 0.078125;
dz=0.078125;

savevtkvector_bin(Vx, Vy, Vz, 'Bshift15000.vtk','B',dx,dy,dz,0,0,0)

end

clear all
close all
yesE=0
if(yesE)
    file='/Volumes/Granata/paraview/tred54/fields_E_cycle15000.vtk'

[V,Ex,Ey,Ez]=read_vtk_3d(file,1);

[Nx,Ny,Nz]=size(Ex);

Vx=circshift(Ex,Nx/2,1);
Vy=circshift(Ey,Nx/2,1);
Vz=circshift(Ez,Nx/2,1);

dx= 0.078125;
dy= 0.078125;
dz=0.078125;

savevtkvector_bin(Vx, Vy, Vz, 'Eshift15000.vtk','E',dx,dy,dz,0,0,0)

end


yesJ=0
if(yesJ)
file='/Volumes/Granata/paraview/tred54/Je_xyz_cycle15000.vtk'

[V,Jx,Jy,Jz]=read_vtk_3d(file,0);

[Nx,Ny,Nz]=size(Jx);

Vx=circshift(Jx,Nx/2,1);
Vy=circshift(Jy,Nx/2,1);
Vz=circshift(Jz,Nx/2,1);

dx= 0.078125;
dy= 0.078125;
dz=0.078125;

savevtkvector_bin(Vx, Vy, Vz, 'Jeshift15000.vtk','Je',dx,dy,dz,0,0,0)

file='/Volumes/Granata/paraview/tred54/Ji_xyz_cycle15000.vtk'

[V,Jx,Jy,Jz]=read_vtk_3d(file,0);

Vx=circshift(Jx,Nx/2,1);
Vy=circshift(Jy,Nx/2,1);
Vz=circshift(Jz,Nx/2,1);

savevtkvector_bin(Vx, Vy, Vz, 'Jishift15000.vtk','Ji',dx,dy,dz,0,0,0)

end