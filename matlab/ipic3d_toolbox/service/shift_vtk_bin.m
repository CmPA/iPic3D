
addpath(genpath('~/matlab4iPic3D/ipic3d_toolbox'));

clear all
close all

directory='/shared/gianni/tred70/'
for cycle=000:1000:24000
ncycle=num2str(cycle)

yesn=1
yesB=1
yesE=1
yesJ=1
yesPe=1
Nsm=0;

if(yesn)
file=[directory 'rho_xyz_cycle' ncycle '.vtk']

[ne,ni]=read_vtk_multiscalar_3d_bin(file,2);

[Nx,Ny,Nz]=size(ne);

ne=circshift(ne,Nx/2);
ni=circshift(ni,Nx/2);
%ne=smooth3Dnew(ne,Nsm);
%ni=smooth3Dnew(ni,Nsm);

dx= 0.078125;
dy= 0.078125;
dz=0.078125;

savevtk_bin(ni, [directory 'Ni_shift' ncycle '.vtk'],'ni',dx,dy,dz,0,0,0)
savevtk_bin(ne, [directory 'Ne_shift' ncycle '.vtk'],'ne',dx,dy,dz,0,0,0)
end


if(yesB)
    file=[directory 'B_xyz_cycle' ncycle '.vtk']

[V,Bx,By,Bz]=read_vtk_3d_bin(file,0);

[Nx,Ny,Nz]=size(Bx);

Vx=circshift(Bx,Nx/2);
Vy=circshift(By,Nx/2);
Vz=circshift(Bz,Nx/2);

dx= 0.078125;
dy= 0.078125;
dz=0.078125;

savevtkvector_bin(Vx, Vy, Vz, [directory 'B_shift' ncycle '.vtk'],'B',dx,dy,dz,0,0,0)

end



if(yesE) 
    file=[directory 'E_xyz_cycle' ncycle '.vtk']

[V,Ex,Ey,Ez]=read_vtk_3d_bin(file,0);


[Nx,Ny,Nz]=size(Ex);

Vx=circshift(Ex,Nx/2); %later versions of matlab require a thrid argumen ,1
Vy=circshift(Ey,Nx/2);
Vz=circshift(Ez,Nx/2);

Vx=smooth3Dnew(Vx,Nsm);
Vy=smooth3Dnew(Vy,Nsm);
Vz=smooth3Dnew(Vz,Nsm);

dx= 0.078125;
dy= 0.078125;
dz=0.078125;

savevtkvector_bin(Vx, Vy, Vz, [directory 'E_shift' ncycle '.vtk'],'E',dx,dy,dz,0,0,0)

end



if(yesJ)
    file=[directory 'Je_xyz_cycle' ncycle '.vtk']

[V,Jx,Jy,Jz]=read_vtk_3d_bin(file,0);

[Nx,Ny,Nz]=size(Jx);

Vx=circshift(Jx,Nx/2);
Vy=circshift(Jy,Nx/2);
Vz=circshift(Jz,Nx/2);

Vx=smooth3Dnew(Vx,Nsm);
Vy=smooth3Dnew(Vy,Nsm);
Vz=smooth3Dnew(Vz,Nsm);

dx= 0.078125;
dy= 0.078125;
dz=0.078125;

savevtkvector_bin(Vx, Vy, Vz, [directory 'Je_shift' ncycle '.vtk'],'Je',dx,dy,dz,0,0,0)


    file=[directory 'Ji_xyz_cycle' ncycle '.vtk']
[V,Jx,Jy,Jz]=read_vtk_3d_bin(file,0);

Vx=circshift(Jx,Nx/2);
Vy=circshift(Jy,Nx/2);
Vz=circshift(Jz,Nx/2);

Vx=smooth3Dnew(Vx,Nsm);
Vy=smooth3Dnew(Vy,Nsm);
Vz=smooth3Dnew(Vz,Nsm);

savevtkvector_bin(Vx, Vy, Vz, [directory 'Ji_shift' ncycle '.vtk'],'Ji',dx,dy,dz,0,0,0)

end


if(yesPe)
    file=[directory 'Pe_xyz_cycle' ncycle '.vtk']

[Pexx,Pexy,Pexz,Peyy,Peyz,Pezz,Pepar,Peper1,Peper2,Peeps,Nx,Ny,Nz] = read_binVTK_pressure(directory,'Pe',cycle);


Pexx=circshift(Pexx,Nx/2);
Pexy=circshift(Pexy,Nx/2);
Pexz=circshift(Pexz,Nx/2);
Peyy=circshift(Peyy,Nx/2);
Peyz=circshift(Peyz,Nx/2);
Pezz=circshift(Pezz,Nx/2);

Pepar=circshift(Pepar,Nx/2);
Peper1=circshift(Peper1,Nx/2);
Peper2=circshift(Peper2,Nx/2);

Nsm=5
Pexx=smooth3Dnew(Pexx,Nsm);
Pexy=smooth3Dnew(Pexy,Nsm);
Pexz=smooth3Dnew(Pexz,Nsm);
Peyy=smooth3Dnew(Peyy,Nsm);
Peyz=smooth3Dnew(Peyz,Nsm);
Pezz=smooth3Dnew(Pezz,Nsm);

Pepar=smooth3Dnew(Pepar,Nsm);
Peper1=smooth3Dnew(Peper1,Nsm);
Peper2=smooth3Dnew(Peper2,Nsm);

dx= 0.078125;
dy= 0.078125;
dz=0.078125;

 file=[directory 'B_xyz_cycle' ncycle '.vtk']

[V,Bx,By,Bz]=read_vtk_3d_bin(file,0);

[Nx,Ny,Nz]=size(Bx);

Bx=circshift(Bx,Nx/2);
By=circshift(By,Nx/2);
Bz=circshift(Bz,Nx/2);

savevtkvector_bin(Pexx, Peyy, Pezz, [directory 'Pe_diag' ncycle '.vtk'],'Pe_diag',dx,dy,dz,0,0,0)
savevtkvector_bin(Pexy, Pexz, Peyz, [directory 'Pe_offdiag' ncycle '.vtk'],'Pe_offdiag',dx,dy,dz,0,0,0)
savevtkvector_bin(Pepar, Peper1, Peper2, [directory 'Pe_mag' ncycle '.vtk'],'Pe_mag',dx,dy,dz,0,0,0)

[Agyro,Agyro_aunai,Nongyro_swisdak]=compute_agyro(Bx, By, Bz, Pexx, Pexy, Pexz, Peyy, Peyz, Pezz ) ;

Nsm=0
Agyro=smooth3Dnew(Agyro,Nsm);
Agyro_aunai=smooth3Dnew(Agyro_aunai,Nsm);
Nongyro_swisdak=smooth3Dnew(Nongyro_swisdak,Nsm);

savevtk_bin(Agyro,[directory 'Agyro_shift' ncycle '.vtk'],'Agyro',dx,dy,dz,0,0,0)
savevtk_bin(Agyro_aunai,[directory 'Agyro_aunai_shift' ncycle '.vtk'],'Agyro_Aunai',dx,dy,dz,0,0,0)
savevtk_bin(Nongyro_swisdak,[directory 'Nongyro_swisdak_shift' ncycle '.vtk'],'Nongyro_swisdak',dx,dy,dz,0,0,0)

        
end

end
