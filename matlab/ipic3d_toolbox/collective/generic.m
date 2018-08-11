addpath(genpath('../ipic3d_toolbox'));

%clearvars -except Ncyc_ini Ncyc_max dir results_dir fraciz Ygsm


dir='~/thr/';
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
Dt=.05;
B0=.0097;
Lx=20;
Ly=60;
Lz=1;
Nx=1280/4;
Ny=3840/4;
Nz=1;
XLEN=20;
YLEN=60;
ZLEN=1;
qom=-256;
mratio=256;
vthi=.045
vthe=.0063
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
