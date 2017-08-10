
addpath(genpath('~/iPic3D/matlab/ipic3d_toolbox'))
%clearvars -except Ncyc_ini Ncyc_max dir results_dir fraciz Ygsm

%dir='/data1/gianni/HRmaha3D3/vtk/'
dir='/data1/gianni/BOW25/'

global Lx Ly Lz Xgsmrange Ygsmrange Zgsmrange dx dy dz XLEN YLEN ZLEN initial_time Nx Ny Nz Dt
Xgsmrange= [-15 -7];
Zgsmrange= [-2.5 2.5];
Ygsmrange= [-1 1];


initial_time=(03*60+48)*60

B0=0.0026

code_E = 240.398;
code_B = 8.01882e-07
%code_B=code_B *1e9; % to convert from Tesla to nT
code_J = 1.20082e-05;
%code_J = code_J*4*pi;
%code_J = code_J*1e9; % to convert to nA/m^2
code_V = 2.99792e+08; %this is c
%code_V=code_V/1e3; %to convert to Km/s
code_T = 2.0468e-12; %this is mp c^2
code_n = 0.25;
code_dp = 5.3140e+04; %ion skin depth for code unit conversion
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
