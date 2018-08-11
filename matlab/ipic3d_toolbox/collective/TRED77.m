addpath(genpath('../ipic3d_toolbox'));

%clearvars -except Ncyc_ini Ncyc_max dir results_dir fraciz Ygsm


dir='/data2/gianni/tred77/';
%dir='/shared/gianni/tred77/'
dir='/Users/Gianni/Dropbox/Science/conferenze/2017/firenze/tred77_15000b/'
global Lx Ly Lz Xgsmrange Ygsmrange Zgsmrange dx dy dz XLEN YLEN ZLEN initial_time Nx Ny Nz Dt



initial_time=0



code_E = 1;
code_B = 1;
code_J = 1;
code_V = 1;
code_T = 1;
code_n = 1;
code_dp = 1; %ion skin depth for code unit conversion
e=1;
mu0 = 4*pi; % to use CGS in code units


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

Xgsmrange= [Lx 0];
Zgsmrange= [0 Ly];
Ygsmrange= [0 Lz];

dx=Lx/Nx;
dy=Ly/Ny;
dz=Lz/Nz;

xc=linspace(0, Lx, Nx);
yc=linspace(0, Ly, Ny);
zc=linspace(0, Lz, Nz);
