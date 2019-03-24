addpath(genpath('../ipic3d_toolbox'));

%clearvars -except Ncyc_ini Ncyc_max dir results_dir fraciz Ygsm


dir='/data2/gianni/turbo/';
global Lx Ly Lz Xgsmrange Ygsmrange Zgsmrange dx dy dz XLEN YLEN ZLEN initial_time Nx Ny Nz Dt



initial_time=0

e=1.6e-19;
mu0 = 4*pi*1e-7;
code_E = 1;
code_B = 1;
%code_B=code_B *1e9; % to convert from Tesla to nT
code_J = 1;
%code_J = code_J*4*pi;
%code_J = code_J*1e9; % to convert to nA/m^2
code_V = 1;
%code_V=code_V/1e3; %to convert to Km/s
code_T = 1;
code_n = 1;
code_dp=1;

filename=[dir 'settings.hdf'];
Dt=1.;
B0=.0097;
Lx=128;
Ly=128;
Lz=0.25;
Nx=512;
Ny=512;
Nz=1;
XLEN=32;
YLEN=32;
ZLEN=1;
qom=-25;
mratio=25;
vthe=.025
vthi=.005
Nprocs=XLEN*YLEN*ZLEN

vthr=vthe;

Xgsmrange= [0 Lx];
Zgsmrange= [0 Ly];
Ygsmrange= [0 Lz];

dx=Lx/Nx;
dy=Ly/Ny;
dz=Lz/Nz;

xc=linspace(0, Lx, Nx);
yc=linspace(0, Ly, Ny);
zc=linspace(0, Lz, Nz);
