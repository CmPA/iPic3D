%%% PLASMA Properties

% Refernce Density in cc
n=1

% Ion Temperature in eV 
T=5e3

% Solar wind Ti/Te

TioTe=5

% B in Tesla
B=20e-9

% DF speed (m/s)
vsw=400e3

% mass ratio code

mrcode=256

% Physics Constants
mu0=4*pi*1e-7;
eps0=8.8542*1.e-12;
c=1/sqrt(mu0*eps0);
k=1.3807e-23;
e= 1.6022e-19;
me=9.1094e-31;
mp=1.6726e-27;

RE= 6371e3; %Earth radius in  meters

% Convert to SI
np=n*1e6; % puts ni in m^-3
ne=np;
Tp=T*1.1604e4; %puts T in K
Te=Tp/TioTe;

% Plasma scales

disp('Protons')
wpp=sqrt(np*e^2/mp/eps0)
dp=c/wpp
vthp=sqrt(k*Tp/mp)
wcp=e*B/mp
rhop=vthp/wcp

disp('Electrons')
wpe=sqrt(ne*e^2/me/eps0)
de=c/wpe
vthe=sqrt(k*Te/me)
wce=e*B/me
rhoe=vthe/wce
lde=sqrt(eps0*k*Te/ne/e^2)

va=B/sqrt(np*mp*mu0)

% Normalisations
%To convert 1 in the code to SI multiply by this
%To convert 1 in SI to code units divide by this
Enorm=c*mp*wpp/e 
Bnorm=mp*wpp/e 
Jnorm=mp*wpp/mu0/e/dp 

