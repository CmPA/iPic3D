function [code_n, code_J, code_V, code_T, code_E, code_B, momentum_corrector] =   code_units()
% Refernce Density in cc
 n_ref = .25;

% Ion Temperature in eV
 Ti_ref = 3.5e3;

% Temperature Ratio
 TioTe = 5;

%B in  Tesla
 B_ref = 20e-9;

%mass ratio code
 mrcode = 25;
% mrcode = 1836;
% boost = 100;

% Physics Constants
 mu0=4*pi*1e-7;
 eps0=8.8542*1.e-12;
 cphys=1/sqrt(mu0*eps0);
 k=1.3807e-23;
 e= 1.6022e-19;
%physical electorns
 me=9.1094e-31;
 mp = me * mrcode;
% physical ions with boost
% mp=1.6726e-27/boost;
% me = mp / mrcode;

% Convert to SI
 np=n_ref*1e6; % puts ni in m^-3
 ne=np;
 Tp=Ti_ref * 1.1604e4; % puts T in K
 Te=Tp/TioTe;

% Plasma scales

%disp('Protons')
 wpp=sqrt(np*e^2/mp/eps0);
 dp=cphys/wpp;
 vthp=sqrt(k*Tp/mp);
 wcp=e*B_ref/mp;
 rhop=vthp/wcp;

%disp('Electrons')
 wpe=sqrt(ne*e^2/me/eps0);
 de=cphys/wpe;
 vthe=sqrt(k*Te/me);
 wce=e*B_ref/me;
 rhoe=vthe/wce;
 lde=sqrt(eps0*k*Te/ne/e*e);

% Normalisations
%To convert 1 in the code to SI multiply by this
%To convert 1 in SI to code units divide by this
code_E = cphys*mp*wpp/e;
code_B = mp*wpp/e;
code_J = mp*wpp/mu0/e/dp;
code_V = cphys;
code_n = n_ref*1e6;
% Assuming the reference species is ions
code_T = cphys *cphys * mp;

momentum_corrector = sqrt(1836.0/mrcode);


return

