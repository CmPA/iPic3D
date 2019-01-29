clear all
load('bimodal_ions_FAC.mat');
close all
addpath(genpath('../../ipic3d_toolbox'));
[Nx,Ny,Nz]=size(fcutoff)
[xg,yg,zg]=ndgrid(0:Nx-1,0:Ny-1,0:Nz-1);
xg=xg/(Nx-1)*560*2-560;
yg=yg/(Ny-1)*560*2-560;
zg=zg/(Nz-1)*560*2-560;
dv= 1/(Nx-1)*560*2;
st=ones(size(fcutoff));
[autX3,autY3,autZ3]=compute_autocorrelation_directional(fcutoff);
[sfX3,sfY3,sfZ3,volX,volY,volZ]=compute_strct_funct_directional(st,fcutoff);
figure(2)
semilogy((1:Nx/2)*dv,sfX3(1:Nx/2))
hold on
semilogy((1:Ny/2)*dv,sfY3(1:Ny/2))
semilogy((1:Nz/2)*dv,sfZ3(1:Nz/2))
legend('X','Y','Z')
xlabel('Km/s')
title('Structure Function')
print('-dpng','str_fun.png')
figure(1)
semilogy((1:Nx)*dv,autX3(1:Nx))
hold on
semilogy((1:Ny)*dv,autY3(1:Ny))
semilogy((1:Nz)*dv,autZ3(1:Nz))
legend('X','Y','Z')
xlabel('Km/s')
title('Autocorrelation')
print('-dpng','autocorrelation.png')