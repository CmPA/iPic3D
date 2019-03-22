cdf_analisys
close all
addpath(genpath('../../ipic3d_toolbox'));
[Nx,Ny,Nz]=size(fcutoff)
st=ones(size(fcutoff));
[autX3,autY3,autZ3]=compute_autocorrelation_directional(fcutoff);
[sfX3,sfY3,sfZ3,volX,volY,volZ]=compute_strct_funct_directional(st,fcutoff);
figure(2)
plot(sfX3)
hold on
plot(sfY3)
plot(sfZ3)
legend('X','Y','Z')
title('Structure Function')
figure(2)
figure(1)
plot(autX3)
hold on
plot(autY3)
plot(autZ3)
legend('X','Y','Z')
title('Autocorrelation')
