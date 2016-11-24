/// here, physical constants used for conversion from code units to IS (at this stage, mostly in the MonteCarlo wrapper)
#include <math.h>


double Realeps0 = 8.8542e-12; // [F/m], permittivity of free space 
double e= 1.6021766208e-19; //[C], elementary charge 
double Realmu0 = 4*M_PI*1e-7;
double mP = 1.6726219e-27; //[Kg]# proton mass  --- atomic mass unit
double k= 1.38064852e-23; //[m2 kg s-2 K-1], Boltzmann constant 
double Realc= 1/sqrt(Realeps0*Realmu0); // [m/s], speed of light   
double R= 8.3144598;// [J K-1 mol-1], gas constant in IS units  
double JtoeV= 6.242e+18; // to convert from Joule to eV 
