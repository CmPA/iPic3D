clear all
close all
load('bimodal_ions_FAC.mat');
[Nx,Ny,Nz]=size(fcutoff);
a=reshape(fcutoff,Nx*Ny*Nz,1);
[pdf,f]=hist(a,1000);
pdf=pdf./sum(pdf);
cdf=fliplr(cumsum(fliplr(pdf)));
loglog(f,cdf)