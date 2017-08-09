addpath(genpath('../ipic3d_toolbox'));

%clearvars -except Ncyc_ini Ncyc_max dir results_dir fraciz Ygsm

%dir='/data1/gianni/HRmaha3D3/vtk/'
dir='/data1/gianni/HRmaha3D3/data3/'

global Lx Ly Lz Xgsmrange Ygsmrange Zgsmrange dx dy dz XLEN YLEN ZLEN initial_time Nx Ny Nz Dt
Xgsmrange= [-45 -15];
Zgsmrange= [-8.7 3.3];
Ygsmrange= [-3 9];

P4 = [-9.8, 3.6, -2.2]

[x_P4 y_P4 z_P4 ix_P4 iy_P4 iz_P4 ipx ipy ipz ip_P4] = gsm2code(P4(1),P4(2),P4(3));

initial_time=(03*60+48)*60

B0=0.0026

code_E = 2060.21;
code_B = 6.87213e-06;
%code_B=code_B *1e9; % to convert from Tesla to nT
code_J = 1.20082e-05;
%code_J = code_J*4*pi;
%code_J = code_J*1e9; % to convert to nA/m^2
code_V = 2.99792e+08;
%code_V=code_V/1e3; %to convert to Km/s
code_T = 1.50326e-10;
code_n = 0.25;
code_dp = 4.5541e+05; %ion skin depth for code unit conversion
e=1.6e-19;
mu0 = 4*pi*1e-7;
%convert to keV
%TeoTi=1/5;
%code_T=code_T/e/1e3/TeoTi;

filename=[dir 'settings.hdf'];
Dt=double(hdf5read(filename,'/collective/Dt'));
B0=double(hdf5read(filename,'collective/Bx0'));
Lx=double(hdf5read(filename,'/collective/Lx'));
Ly=double(hdf5read(filename,'/collective/Ly'));
Lz=double(hdf5read(filename,'/collective/Lz'));
Nx=double(hdf5read(filename,'/collective/Nxc'));
Ny=double(hdf5read(filename,'/collective/Nyc'));
Nz=double(hdf5read(filename,'/collective/Nzc'));
XLEN=double(hdf5read(filename,'/topology/XLEN'));
YLEN=double(hdf5read(filename,'/topology/YLEN'));
ZLEN=double(hdf5read(filename,'/topology/ZLEN'));
qom=double(hdf5read(filename,'/collective/species_0/qom'))
mratio=abs(hdf5read(filename,'/collective/species_0/qom'));
vthi=double(hdf5read(filename,'/collective/species_1/uth'))
vthe=double(hdf5read(filename,'/collective/species_0/uth'))
Nprocs=hdf5read(filename,'/topology/Nprocs')

vthr=vthe;

dx=Lx/Nx;
dy=Ly/Ny;
dz=Lz/Nz;

xc=linspace(0, Lx, Nx);
yc=linspace(0, Ly, Ny);
zc=linspace(0, Lz, Nz);
