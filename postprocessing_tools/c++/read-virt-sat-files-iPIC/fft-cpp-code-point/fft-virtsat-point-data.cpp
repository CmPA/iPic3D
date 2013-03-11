/***************************************************************************

    fft-virtsat-point-data.cpp  -  FFT of 1D real data extracted with process-virtual-sat-point.cpp

    Author:  A.E.Vapirev, KU Leuven, Afdeling Plasma-astrofysica
             2011 Oct 7

    Output               : The absolute value of the FFT is being written out and the x-axis is frequency
                           in radians/wci which are the frequencies correspoinding to the FFT elements.
                           Then this frequency is scaled - see the Note below.

    Change log           :

    If using FFTW3, compile with "c++ fft-virtsat-point.cpp -I/usr/local/include/ -L/usr/lib/ -lfftw3 -lm"
                         or with "c++ fft-virtsat-point.cpp -lfftw3 -lm"
                         or if you install FFTW3 from source adjust the 'include' and 'lib' locations accordingly.
                         See here for examples http://opencv-code.com/Using_the_Fast_Fourier_Transform_Library_FFTW

    Note                 : The following quantities are needed for scaling the FFT output x-axis

                           n0=mean(rhoi-rhoe)/2;
                           b0=sqrt(mean(bx.^2+by.^2+bz.^2));
                           wci=b0;
                           wpi=1*sqrt(n0);
                           wce=wci*mratio;
                           wlh=1/sqrt(1/wce/wci+1/wpi^2);

                           The FFT contains information between 0 and fs, where fs is the sample frequency.
                           The sampling frequency is the frequency rate between the adjacent samples.
                           It is calculated as fs=1/delta_T. In the case of iPIC3D delta_T=wci*dt
                           So if the frequency axis is plotted as (0:(M-1))/M it is actually f/fs what is plotted over X-axis.
                           If we multiply by 2Pi/delta_T, e.g. [(0:(M-1))/M]*[2Pi/delta_T] we get it in radians/wci.
                           If we multiply by 2Pi/dt, e.g. [(0:(M-1))/M]*[2Pi/dt] we get it in radians.


 ***************************************************************************/                      

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <complex>
#include "fftw3.h"

int main(int argc, char* argv[]) {

    int mass_ratio;
    int istart_record;
    int istop_record;
    float B0x;
    float dt;

    sscanf(argv[1],"%d",&mass_ratio);
    sscanf(argv[2],"%d",&istart_record);
    sscanf(argv[3],"%d",&istop_record);
    sscanf(argv[4],"%f",&B0x);
    sscanf(argv[5],"%f",&dt);

    int mratio = abs(mass_ratio);
    int istart = istart_record;
    int istop  = istop_record;

    // get the number of data points from 'BestVirtualSatelliteInfoFoundPoint.txt'
    char VirtualSatelliteFilenameFound[] = "";  // dummy
    int satellite_number_best_match;            // dummy
    int number_of_data_rows;
    char info_fname[] = "BestVirtualSatelliteInfoFoundPoint.txt";
    FILE *satellite_info_file;
    satellite_info_file = fopen(info_fname,"r");
    fscanf(satellite_info_file,"%s %i %i",VirtualSatelliteFilenameFound,&satellite_number_best_match,&number_of_data_rows);
    fclose(satellite_info_file);

    int N=number_of_data_rows;   // number of input data points to read
    int NFFT=istop-istart;       // number of data elements (data window) on which we perform the FFTs
    int M=NFFT/2+1;              // number of output FFT data points

    // quantities to read from the file
    float omega_ci_x_t;

    float Bx[number_of_data_rows];
    float By[number_of_data_rows];
    float Bz[number_of_data_rows];
    float Ex[number_of_data_rows];
    float Ey[number_of_data_rows];
    float Ez[number_of_data_rows];
    float Jx_tot_e[number_of_data_rows];
    float Jy_tot_e[number_of_data_rows];
    float Jz_tot_e[number_of_data_rows];
    float Jx_tot_i[number_of_data_rows];
    float Jy_tot_i[number_of_data_rows];
    float Jz_tot_i[number_of_data_rows];
    float rho_tot_e[number_of_data_rows];
    float rho_tot_i[number_of_data_rows];
    float Vx_tot_e[number_of_data_rows];
    float Vy_tot_e[number_of_data_rows];
    float Vz_tot_e[number_of_data_rows];
    float Vx_tot_i[number_of_data_rows];
    float Vy_tot_i[number_of_data_rows];
    float Vz_tot_i[number_of_data_rows];

    float E[number_of_data_rows];
    float B[number_of_data_rows];
    float V_tot_e[number_of_data_rows];
    float V_tot_i[number_of_data_rows];

    float dBdt[number_of_data_rows];
    float dEdt[number_of_data_rows];

    // define input for FFT solver
    double * FFTin_Bx  = (double*) malloc(sizeof(double) * NFFT);
    double * FFTin_By  = (double*) malloc(sizeof(double) * NFFT);
    double * FFTin_Bz  = (double*) malloc(sizeof(double) * NFFT);
    double * FFTin_Ex  = (double*) malloc(sizeof(double) * NFFT);
    double * FFTin_Ey  = (double*) malloc(sizeof(double) * NFFT);
    double * FFTin_Ez  = (double*) malloc(sizeof(double) * NFFT);
    double * FFTin_Jx_tot_e  = (double*) malloc(sizeof(double) * NFFT);
    double * FFTin_Jy_tot_e  = (double*) malloc(sizeof(double) * NFFT);
    double * FFTin_Jz_tot_e  = (double*) malloc(sizeof(double) * NFFT);
    double * FFTin_Jx_tot_i  = (double*) malloc(sizeof(double) * NFFT);
    double * FFTin_Jy_tot_i  = (double*) malloc(sizeof(double) * NFFT);
    double * FFTin_Jz_tot_i  = (double*) malloc(sizeof(double) * NFFT);
    double * FFTin_rho_tot_e  = (double*) malloc(sizeof(double) * NFFT);
    double * FFTin_rho_tot_i  = (double*) malloc(sizeof(double) * NFFT);
    double * FFTin_Vx_tot_e  = (double*) malloc(sizeof(double) * NFFT);
    double * FFTin_Vy_tot_e  = (double*) malloc(sizeof(double) * NFFT);
    double * FFTin_Vz_tot_e  = (double*) malloc(sizeof(double) * NFFT);
    double * FFTin_Vx_tot_i  = (double*) malloc(sizeof(double) * NFFT);
    double * FFTin_Vy_tot_i  = (double*) malloc(sizeof(double) * NFFT);
    double * FFTin_Vz_tot_i  = (double*) malloc(sizeof(double) * NFFT);
    double * FFTin_B  = (double*) malloc(sizeof(double) * NFFT);
    double * FFTin_E  = (double*) malloc(sizeof(double) * NFFT);
    double * FFTin_V_tot_e  = (double*) malloc(sizeof(double) * NFFT);
    double * FFTin_V_tot_i  = (double*) malloc(sizeof(double) * NFFT);
    // define output for FFT solver
    fftw_complex * FFTout_Bx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    fftw_complex * FFTout_By = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    fftw_complex * FFTout_Bz = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    fftw_complex * FFTout_Ex = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    fftw_complex * FFTout_Ey = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    fftw_complex * FFTout_Ez = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    fftw_complex * FFTout_Jx_tot_e = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    fftw_complex * FFTout_Jy_tot_e = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    fftw_complex * FFTout_Jz_tot_e = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    fftw_complex * FFTout_Jx_tot_i = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    fftw_complex * FFTout_Jy_tot_i = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    fftw_complex * FFTout_Jz_tot_i = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    fftw_complex * FFTout_rho_tot_e = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    fftw_complex * FFTout_rho_tot_i = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    fftw_complex * FFTout_Vx_tot_e = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    fftw_complex * FFTout_Vy_tot_e = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    fftw_complex * FFTout_Vz_tot_e = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    fftw_complex * FFTout_Vx_tot_i = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    fftw_complex * FFTout_Vy_tot_i = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    fftw_complex * FFTout_Vz_tot_i = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    fftw_complex * FFTout_B  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    fftw_complex * FFTout_E  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    fftw_complex * FFTout_V_tot_e  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    fftw_complex * FFTout_V_tot_i  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
    // define plan for FFT solver
    fftw_plan plan_Bx = fftw_plan_dft_r2c_1d(NFFT, FFTin_Bx, FFTout_Bx, FFTW_ESTIMATE);
    fftw_plan plan_By = fftw_plan_dft_r2c_1d(NFFT, FFTin_By, FFTout_By, FFTW_ESTIMATE);
    fftw_plan plan_Bz = fftw_plan_dft_r2c_1d(NFFT, FFTin_Bz, FFTout_Bz, FFTW_ESTIMATE);
    fftw_plan plan_Ex = fftw_plan_dft_r2c_1d(NFFT, FFTin_Ex, FFTout_Ex, FFTW_ESTIMATE);
    fftw_plan plan_Ey = fftw_plan_dft_r2c_1d(NFFT, FFTin_Ey, FFTout_Ey, FFTW_ESTIMATE);
    fftw_plan plan_Ez = fftw_plan_dft_r2c_1d(NFFT, FFTin_Ez, FFTout_Ez, FFTW_ESTIMATE);
    fftw_plan plan_Jx_tot_e = fftw_plan_dft_r2c_1d(NFFT, FFTin_Jx_tot_e, FFTout_Jx_tot_e, FFTW_ESTIMATE);
    fftw_plan plan_Jy_tot_e = fftw_plan_dft_r2c_1d(NFFT, FFTin_Jy_tot_e, FFTout_Jy_tot_e, FFTW_ESTIMATE);
    fftw_plan plan_Jz_tot_e = fftw_plan_dft_r2c_1d(NFFT, FFTin_Jz_tot_e, FFTout_Jz_tot_e, FFTW_ESTIMATE);
    fftw_plan plan_Jx_tot_i = fftw_plan_dft_r2c_1d(NFFT, FFTin_Jx_tot_i, FFTout_Jx_tot_i, FFTW_ESTIMATE);
    fftw_plan plan_Jy_tot_i = fftw_plan_dft_r2c_1d(NFFT, FFTin_Jy_tot_i, FFTout_Jy_tot_i, FFTW_ESTIMATE);
    fftw_plan plan_Jz_tot_i = fftw_plan_dft_r2c_1d(NFFT, FFTin_Jz_tot_i, FFTout_Jz_tot_i, FFTW_ESTIMATE);
    fftw_plan plan_rho_tot_e = fftw_plan_dft_r2c_1d(NFFT, FFTin_rho_tot_e, FFTout_rho_tot_e, FFTW_ESTIMATE);
    fftw_plan plan_rho_tot_i = fftw_plan_dft_r2c_1d(NFFT, FFTin_rho_tot_i, FFTout_rho_tot_i, FFTW_ESTIMATE);
    fftw_plan plan_Vx_tot_e = fftw_plan_dft_r2c_1d(NFFT, FFTin_Vx_tot_e, FFTout_Vx_tot_e, FFTW_ESTIMATE);
    fftw_plan plan_Vy_tot_e = fftw_plan_dft_r2c_1d(NFFT, FFTin_Vy_tot_e, FFTout_Vy_tot_e, FFTW_ESTIMATE);
    fftw_plan plan_Vz_tot_e = fftw_plan_dft_r2c_1d(NFFT, FFTin_Vz_tot_e, FFTout_Vz_tot_e, FFTW_ESTIMATE);
    fftw_plan plan_Vx_tot_i = fftw_plan_dft_r2c_1d(NFFT, FFTin_Vx_tot_i, FFTout_Vx_tot_i, FFTW_ESTIMATE);
    fftw_plan plan_Vy_tot_i = fftw_plan_dft_r2c_1d(NFFT, FFTin_Vy_tot_i, FFTout_Vy_tot_i, FFTW_ESTIMATE);
    fftw_plan plan_Vz_tot_i = fftw_plan_dft_r2c_1d(NFFT, FFTin_Vz_tot_i, FFTout_Vz_tot_i, FFTW_ESTIMATE);
    fftw_plan plan_B = fftw_plan_dft_r2c_1d(NFFT, FFTin_B, FFTout_B, FFTW_ESTIMATE);
    fftw_plan plan_E = fftw_plan_dft_r2c_1d(NFFT, FFTin_E, FFTout_E, FFTW_ESTIMATE);
    fftw_plan plan_V_tot_e = fftw_plan_dft_r2c_1d(NFFT, FFTin_V_tot_e, FFTout_V_tot_e, FFTW_ESTIMATE);
    fftw_plan plan_V_tot_i = fftw_plan_dft_r2c_1d(NFFT, FFTin_V_tot_i, FFTout_V_tot_i, FFTW_ESTIMATE);

    float n0 = 0.0;
    float b0 = 0.0;
    float wci;
    float wpi;
    float wce;
    float wlh;
    float wsm;

    // read from input data file
    char data_fname[] = "sat_output_point.txt";
    FILE *input_file;
    input_file = fopen(data_fname,"r");

    // skip (void read) first line (header) in file
    int idummy;
    do idummy = fgetc(input_file); while (idummy != '\n');

    // skip first istart-1 data records in file
    for (int k=0;k<istart;k++){
        do idummy = fgetc(input_file); while (idummy != '\n');
    }

    for (int i=istart;i<istop+1;i++){                                // loop time records
        // read quantities
        fscanf(input_file,"%f ",&omega_ci_x_t);
        fscanf(input_file,"%E %E %E %E %E %E %E %E %E %E %E %E %E %E "
                         ,&Bx[i],&By[i],&Bz[i]
                         ,&Ex[i],&Ey[i],&Ez[i]
                         ,&Jx_tot_e[i],&Jy_tot_e[i],&Jz_tot_e[i]
                         ,&Jx_tot_i[i],&Jy_tot_i[i],&Jz_tot_i[i]
                         ,&rho_tot_e[i],&rho_tot_i[i]);
        fscanf(input_file,"%E %E %E %E %E %E "
                         ,&Vx_tot_e[i],&Vy_tot_e[i],&Vz_tot_e[i]
                         ,&Vx_tot_i[i],&Vy_tot_i[i],&Vz_tot_i[i]);
        fscanf(input_file,"%E %E "
                         ,&B[i],&E[i]);
        fscanf(input_file,"%E %E "
                         ,&V_tot_e[i],&V_tot_i[i]);
        // set real part for FFTs
        FFTin_Bx[i-istart] = Bx[i];
        FFTin_By[i-istart] = By[i];
        FFTin_Bz[i-istart] = Bz[i];
        FFTin_Ex[i-istart] = Ex[i];
        FFTin_Ey[i-istart] = Ey[i];
        FFTin_Ez[i-istart] = Ez[i];
        FFTin_Jx_tot_e[i-istart] = Jx_tot_e[i];
        FFTin_Jy_tot_e[i-istart] = Jy_tot_e[i];
        FFTin_Jz_tot_e[i-istart] = Jz_tot_e[i];
        FFTin_Jx_tot_i[i-istart] = Jx_tot_i[i];
        FFTin_Jy_tot_i[i-istart] = Jy_tot_i[i];
        FFTin_Jz_tot_i[i-istart] = Jz_tot_i[i];
        FFTin_rho_tot_e[i-istart] = rho_tot_e[i];
        FFTin_rho_tot_i[i-istart] = rho_tot_i[i];
        FFTin_Vx_tot_e[i-istart] = Vx_tot_e[i];
        FFTin_Vy_tot_e[i-istart] = Vy_tot_e[i];
        FFTin_Vz_tot_e[i-istart] = Vz_tot_e[i];
        FFTin_Vx_tot_i[i-istart] = Vx_tot_i[i];
        FFTin_Vy_tot_i[i-istart] = Vy_tot_i[i];
        FFTin_Vz_tot_i[i-istart] = Vz_tot_i[i];
        FFTin_B[i-istart] = B[i];
        FFTin_E[i-istart] = E[i];
        FFTin_V_tot_e[i-istart] = V_tot_e[i];
        FFTin_V_tot_i[i-istart] = V_tot_i[i];

        n0 = n0 + rho_tot_i[i]-rho_tot_e[i];                            // '-' for rho_electron because of their negative charge
        b0 = b0 + sqrt( Bx[i]*Bx[i] + By[i]*By[i] + Bz[i]*Bz[i] );
    }

    // close input file
    fclose(input_file);

    /* compute the reference frequencies
    
       wpi - plasma frequency            : [cgs] units?
       wci - ion cyclotron frequency     : wci=charge*b0/m_ion. In iPIC3D charge/m_ion=1 --> wci=b0
       wce - electron cyclotron frequency: wce=wci*mass_ratio
       wlh - lower hybrid frequency      :
       wsm - square mass ratio frequency :

       there are three ways to compute them:
       1. take only first time record for n0 and b0 for the selected interval
       2. take first and last record for n0 and b0 for the selected interval
       3. take averaged values of n0 and b0 over the selected interval

       Seems like taking the first record only fits best the results from the matlab FFT code,
       although taking averaged vaules is the logical way to it. */

/*    // compute reference frequencies
    n0 = n0*4.0*3.1415927/2.0/(istop-istart);                             // compute average and *4Pi in order to convert to real charge density *** something is wrong - gives negative close to zero: - matlab code is n0=mean(rhoi-rhoe)/2;
//    n0 = ( rho_tot_i[istart] - rho_tot_e[istart] )*4.0*3.1415927/2.0;     // for first time record only
    b0 = b0/(istop-istart);                                               // compute the average
//    b0=B0x;
    wci = b0;                                                             // check what frequency is equal to below - matlab code is wci=b0=sqrt(mean(bx.^2+by.^2+bz.^2));
    wpi = sqrt(fabs(n0));                                                 // use with the averaged n0 above
//    wpi = sqrt(n0);
    wce = wci*mratio;
    wlh = 1.0/sqrt(1.0/wce/wci+1.0/wpi/wpi);
    wsm = sqrt(mratio)*wpi;*/ 

    // compute reference frequencies for first time record only
    n0 = ( rho_tot_i[istart] - rho_tot_e[istart] )*4.0*3.1415927/2.0;     // '-' for rho_electron because of their negative charge
    b0 = sqrt( Bx[istart]*Bx[istart] + By[istart]*By[istart] + Bz[istart]*Bz[istart] );
    wci = B0x; //B0x;
    wpi = sqrt(n0);
    wce = wci*mratio;
    wlh = 1.0/sqrt(1.0/wce/wci+1.0/wpi/wpi);
    wsm = sqrt(mratio)*wpi;

//  debugging
    std::cout << "b0 = " << b0 << " " << "n0 = " << n0 << " " << std::endl;
    std::cout << "wci = " << wci << " " << std::endl;
    std::cout << "wpi = " << wpi << " " << std::endl;
    std::cout << "wce = " << wce << " " << std::endl;
    std::cout << "wlh = " << wlh << " " << std::endl; 
    std::cout << "wsm = " << wsm << " " << std::endl;

    // write five reference frequencies to compare later with the FFT output
    char refreq_fname[] = "reference_frequencies_fft_sat_output_point.txt";
    FILE *refreq_file;
    refreq_file = fopen(refreq_fname,"w");
    fprintf(refreq_file,"%f %f %f %f %f ",wci,wpi,wce,wlh,wsm);
    fclose(refreq_file);

    // perform the FFT
    fftw_execute(plan_Bx);
    fftw_execute(plan_By);
    fftw_execute(plan_Bz);
    fftw_execute(plan_Ex);
    fftw_execute(plan_Ey);
    fftw_execute(plan_Ez);
    fftw_execute(plan_Jx_tot_e);
    fftw_execute(plan_Jy_tot_e);
    fftw_execute(plan_Jz_tot_e);
    fftw_execute(plan_Jx_tot_i);
    fftw_execute(plan_Jy_tot_i);
    fftw_execute(plan_Jz_tot_i);
    fftw_execute(plan_rho_tot_e);
    fftw_execute(plan_rho_tot_i);
    fftw_execute(plan_Vx_tot_e);
    fftw_execute(plan_Vy_tot_e);
    fftw_execute(plan_Vz_tot_e);
    fftw_execute(plan_Vx_tot_i);
    fftw_execute(plan_Vy_tot_i);
    fftw_execute(plan_Vz_tot_i);
    fftw_execute(plan_B);
    fftw_execute(plan_E);
    fftw_execute(plan_V_tot_e);
    fftw_execute(plan_V_tot_i);

    // write out the FFTs

    // open the output file
    char output_filename[] = "fft_sat_output_point.txt";
    FILE *output_file;
    output_file = fopen(output_filename,"w");

    // write out notation at the beginning
    fprintf(output_file,"%s","FFTs of  ");
    fprintf(output_file,"%s","Bx            By            Bz            Ex            Ey            Ez");
    fprintf(output_file,"%s","     Jx_tot_e      Jy_tot_e      Jz_tot_e      Jx_tot_i");
    fprintf(output_file,"%s","      Jy_tot_i      Jz_tot_i    rho_tot_e    rho_tot_i      ");
    fprintf(output_file,"%s","Vx_tot_e     Vy_tot_e     Vz_tot_e     Vx_tot_i     Vy_tot_i     Vz_tot_i");
    fprintf(output_file,"%s","          B            E            ");
    fprintf(output_file,"%s","   V_tot_e     V_tot_i");
    fprintf(output_file,"%s","dBdt         dEdt          dBxdt        dBydt        dBzdt        dExdt        dEydt        dEzdt      ");
    fprintf(output_file,"\n","");

    // define abs values of the FFT output
    float absFFTout_Bx;
    float absFFTout_By;
    float absFFTout_Bz;
    float absFFTout_Ex;
    float absFFTout_Ey;
    float absFFTout_Ez;
    float absFFTout_Jx_tot_e;
    float absFFTout_Jy_tot_e;
    float absFFTout_Jz_tot_e;
    float absFFTout_Jx_tot_i;
    float absFFTout_Jy_tot_i;
    float absFFTout_Jz_tot_i;
    float absFFTout_rho_tot_e;
    float absFFTout_rho_tot_i;
    float absFFTout_Vx_tot_e;
    float absFFTout_Vy_tot_e;
    float absFFTout_Vz_tot_e;
    float absFFTout_Vx_tot_i;
    float absFFTout_Vy_tot_i;
    float absFFTout_Vz_tot_i;
    float absFFTout_B;
    float absFFTout_E;
    float absFFTout_V_tot_e;
    float absFFTout_V_tot_i;

/*    // find abs values of the FFT output
    for (int i=0;i<M;i++){                                // loop time records
        absFFTout_Bx=sqrt(FFTout_Bx[i][0]*FFTout_Bx[i][0]+FFTout_Bx[i][1]*FFTout_Bx[i][1]);
        absFFTout_By=sqrt(FFTout_By[i][0]*FFTout_By[i][0]+FFTout_By[i][1]*FFTout_By[i][1]);
        absFFTout_Bz=sqrt(FFTout_Bz[i][0]*FFTout_Bz[i][0]+FFTout_Bz[i][1]*FFTout_Bz[i][1]);
        absFFTout_Ex=sqrt(FFTout_Ex[i][0]*FFTout_Ex[i][0]+FFTout_Ex[i][1]*FFTout_Ex[i][1]);
        absFFTout_Ey=sqrt(FFTout_Ey[i][0]*FFTout_Ey[i][0]+FFTout_Ey[i][1]*FFTout_Ey[i][1]);
        absFFTout_Ez=sqrt(FFTout_Ez[i][0]*FFTout_Ez[i][0]+FFTout_Ez[i][1]*FFTout_Ez[i][1]);
        absFFTout_Jx_tot_e=sqrt(FFTout_Jx_tot_e[i][0]*FFTout_Jx_tot_e[i][0]+FFTout_Jx_tot_e[i][1]*FFTout_Jx_tot_e[i][1]);
        absFFTout_Jy_tot_e=sqrt(FFTout_Jy_tot_e[i][0]*FFTout_Jy_tot_e[i][0]+FFTout_Jy_tot_e[i][1]*FFTout_Jy_tot_e[i][1]);
        absFFTout_Jz_tot_e=sqrt(FFTout_Jz_tot_e[i][0]*FFTout_Jz_tot_e[i][0]+FFTout_Jz_tot_e[i][1]*FFTout_Jz_tot_e[i][1]);
        absFFTout_Jx_tot_i=sqrt(FFTout_Jx_tot_i[i][0]*FFTout_Jx_tot_i[i][0]+FFTout_Jx_tot_i[i][1]*FFTout_Jx_tot_i[i][1]);
        absFFTout_Jy_tot_i=sqrt(FFTout_Jy_tot_i[i][0]*FFTout_Jy_tot_i[i][0]+FFTout_Jy_tot_i[i][1]*FFTout_Jy_tot_i[i][1]);
        absFFTout_Jz_tot_i=sqrt(FFTout_Jz_tot_i[i][0]*FFTout_Jz_tot_i[i][0]+FFTout_Jz_tot_i[i][1]*FFTout_Jz_tot_i[i][1]);
        absFFTout_rho_tot_e=sqrt(FFTout_rho_tot_e[i][0]*FFTout_rho_tot_e[i][0]+FFTout_rho_tot_e[i][1]*FFTout_rho_tot_e[i][1]);
        absFFTout_rho_tot_i=sqrt(FFTout_rho_tot_i[i][0]*FFTout_rho_tot_i[i][0]+FFTout_rho_tot_i[i][1]*FFTout_rho_tot_i[i][1]);
        absFFTout_Vx_tot_e=sqrt(FFTout_Vx_tot_e[i][0]*FFTout_Vx_tot_e[i][0]+FFTout_Vx_tot_e[i][1]*FFTout_Vx_tot_e[i][1]);
        absFFTout_Vy_tot_e=sqrt(FFTout_Vy_tot_e[i][0]*FFTout_Vy_tot_e[i][0]+FFTout_Vy_tot_e[i][1]*FFTout_Vy_tot_e[i][1]);
        absFFTout_Vz_tot_e=sqrt(FFTout_Vz_tot_e[i][0]*FFTout_Vz_tot_e[i][0]+FFTout_Vz_tot_e[i][1]*FFTout_Vz_tot_e[i][1]);
        absFFTout_Vx_tot_i=sqrt(FFTout_Vx_tot_i[i][0]*FFTout_Vx_tot_i[i][0]+FFTout_Vx_tot_i[i][1]*FFTout_Vx_tot_i[i][1]);
        absFFTout_Vy_tot_i=sqrt(FFTout_Vy_tot_i[i][0]*FFTout_Vy_tot_i[i][0]+FFTout_Vy_tot_i[i][1]*FFTout_Vy_tot_i[i][1]);
        absFFTout_Vz_tot_i=sqrt(FFTout_Vz_tot_i[i][0]*FFTout_Vz_tot_i[i][0]+FFTout_Vz_tot_i[i][1]*FFTout_Vz_tot_i[i][1]);
        absFFTout_B=sqrt(FFTout_B[i][0]*FFTout_B[i][0]+FFTout_B[i][1]*FFTout_B[i][1]);
        absFFTout_E=sqrt(FFTout_E[i][0]*FFTout_E[i][0]+FFTout_E[i][1]*FFTout_E[i][1]);
        absFFTout_V_tot_e=sqrt(FFTout_V_tot_e[i][0]*FFTout_V_tot_e[i][0]+FFTout_V_tot_e[i][1]*FFTout_V_tot_e[i][1]);
        absFFTout_V_tot_i=sqrt(FFTout_V_tot_i[i][0]*FFTout_V_tot_i[i][0]+FFTout_V_tot_i[i][1]*FFTout_V_tot_i[i][1]);
    }*/

    // the x-axis is frequency - see the note in the beginning
    float frequency_radians,frequency_over_Fs;

    // write out the data
    for (int i=0;i<M;i++){                                    // loop time records - time units are wci*t (cyclotron periods)
        frequency_over_Fs=1.0*i/M;
        frequency_radians=frequency_over_Fs;                  // the x-axis is f/Fs - see the note in the beginning
        frequency_radians=2.0*3.14159*frequency_over_Fs/dt;   // the x-axis is [rad] - see the note in the beginning

        // find abs of FFT of the FFT output
        absFFTout_Bx=sqrt(FFTout_Bx[i][0]*FFTout_Bx[i][0]+FFTout_Bx[i][1]*FFTout_Bx[i][1]);
        absFFTout_By=sqrt(FFTout_By[i][0]*FFTout_By[i][0]+FFTout_By[i][1]*FFTout_By[i][1]);
        absFFTout_Bz=sqrt(FFTout_Bz[i][0]*FFTout_Bz[i][0]+FFTout_Bz[i][1]*FFTout_Bz[i][1]);
        absFFTout_Ex=sqrt(FFTout_Ex[i][0]*FFTout_Ex[i][0]+FFTout_Ex[i][1]*FFTout_Ex[i][1]);
        absFFTout_Ey=sqrt(FFTout_Ey[i][0]*FFTout_Ey[i][0]+FFTout_Ey[i][1]*FFTout_Ey[i][1]);
        absFFTout_Ez=sqrt(FFTout_Ez[i][0]*FFTout_Ez[i][0]+FFTout_Ez[i][1]*FFTout_Ez[i][1]);
        absFFTout_Jx_tot_e=sqrt(FFTout_Jx_tot_e[i][0]*FFTout_Jx_tot_e[i][0]+FFTout_Jx_tot_e[i][1]*FFTout_Jx_tot_e[i][1]);
        absFFTout_Jy_tot_e=sqrt(FFTout_Jy_tot_e[i][0]*FFTout_Jy_tot_e[i][0]+FFTout_Jy_tot_e[i][1]*FFTout_Jy_tot_e[i][1]);
        absFFTout_Jz_tot_e=sqrt(FFTout_Jz_tot_e[i][0]*FFTout_Jz_tot_e[i][0]+FFTout_Jz_tot_e[i][1]*FFTout_Jz_tot_e[i][1]);
        absFFTout_Jx_tot_i=sqrt(FFTout_Jx_tot_i[i][0]*FFTout_Jx_tot_i[i][0]+FFTout_Jx_tot_i[i][1]*FFTout_Jx_tot_i[i][1]);
        absFFTout_Jy_tot_i=sqrt(FFTout_Jy_tot_i[i][0]*FFTout_Jy_tot_i[i][0]+FFTout_Jy_tot_i[i][1]*FFTout_Jy_tot_i[i][1]);
        absFFTout_Jz_tot_i=sqrt(FFTout_Jz_tot_i[i][0]*FFTout_Jz_tot_i[i][0]+FFTout_Jz_tot_i[i][1]*FFTout_Jz_tot_i[i][1]);
        absFFTout_rho_tot_e=sqrt(FFTout_rho_tot_e[i][0]*FFTout_rho_tot_e[i][0]+FFTout_rho_tot_e[i][1]*FFTout_rho_tot_e[i][1]);
        absFFTout_rho_tot_i=sqrt(FFTout_rho_tot_i[i][0]*FFTout_rho_tot_i[i][0]+FFTout_rho_tot_i[i][1]*FFTout_rho_tot_i[i][1]);
        absFFTout_Vx_tot_e=sqrt(FFTout_Vx_tot_e[i][0]*FFTout_Vx_tot_e[i][0]+FFTout_Vx_tot_e[i][1]*FFTout_Vx_tot_e[i][1]);
        absFFTout_Vy_tot_e=sqrt(FFTout_Vy_tot_e[i][0]*FFTout_Vy_tot_e[i][0]+FFTout_Vy_tot_e[i][1]*FFTout_Vy_tot_e[i][1]);
        absFFTout_Vz_tot_e=sqrt(FFTout_Vz_tot_e[i][0]*FFTout_Vz_tot_e[i][0]+FFTout_Vz_tot_e[i][1]*FFTout_Vz_tot_e[i][1]);
        absFFTout_Vx_tot_i=sqrt(FFTout_Vx_tot_i[i][0]*FFTout_Vx_tot_i[i][0]+FFTout_Vx_tot_i[i][1]*FFTout_Vx_tot_i[i][1]);
        absFFTout_Vy_tot_i=sqrt(FFTout_Vy_tot_i[i][0]*FFTout_Vy_tot_i[i][0]+FFTout_Vy_tot_i[i][1]*FFTout_Vy_tot_i[i][1]);
        absFFTout_Vz_tot_i=sqrt(FFTout_Vz_tot_i[i][0]*FFTout_Vz_tot_i[i][0]+FFTout_Vz_tot_i[i][1]*FFTout_Vz_tot_i[i][1]);
        absFFTout_B=sqrt(FFTout_B[i][0]*FFTout_B[i][0]+FFTout_B[i][1]*FFTout_B[i][1]);
        absFFTout_E=sqrt(FFTout_E[i][0]*FFTout_E[i][0]+FFTout_E[i][1]*FFTout_E[i][1]);
        absFFTout_V_tot_e=sqrt(FFTout_V_tot_e[i][0]*FFTout_V_tot_e[i][0]+FFTout_V_tot_e[i][1]*FFTout_V_tot_e[i][1]);
        absFFTout_V_tot_i=sqrt(FFTout_V_tot_i[i][0]*FFTout_V_tot_i[i][0]+FFTout_V_tot_i[i][1]*FFTout_V_tot_i[i][1]);

        fprintf(output_file,"%f ",frequency_radians);
        fprintf(output_file,"%E %E %E %E %E %E %E %E %E %E %E %E %E %E "
                           ,absFFTout_Bx,absFFTout_By,absFFTout_Bz
                           ,absFFTout_Ex,absFFTout_Ey,absFFTout_Ez
                           ,absFFTout_Jx_tot_e,absFFTout_Jy_tot_e,absFFTout_Jz_tot_e
                           ,absFFTout_Jx_tot_i,absFFTout_Jy_tot_i,absFFTout_Jz_tot_i
                           ,absFFTout_rho_tot_e,absFFTout_rho_tot_i);
        fprintf(output_file,"%E %E %E %E %E %E "
                           ,absFFTout_Vx_tot_e,absFFTout_Vy_tot_e,absFFTout_Vz_tot_e
                           ,absFFTout_Vx_tot_i,absFFTout_Vy_tot_i,absFFTout_Vz_tot_i);
        fprintf(output_file,"%E %E "
                           ,absFFTout_B,absFFTout_E);
        fprintf(output_file,"%E %E "
                           ,absFFTout_V_tot_e,absFFTout_V_tot_i);
        fprintf(output_file,"\n","");
    }

    // close the output file
    fclose(output_file);

    // clear the memory from FFT

    fftw_free(FFTout_Bx);
    fftw_free(FFTout_By);
    fftw_free(FFTout_Bz);
    fftw_free(FFTout_Ex);
    fftw_free(FFTout_Ey);
    fftw_free(FFTout_Ez);
    fftw_free(FFTout_Jx_tot_e);
    fftw_free(FFTout_Jy_tot_e);
    fftw_free(FFTout_Jz_tot_e);
    fftw_free(FFTout_Jx_tot_i);
    fftw_free(FFTout_Jy_tot_i);
    fftw_free(FFTout_Jz_tot_i);
    fftw_free(FFTout_rho_tot_e);
    fftw_free(FFTout_rho_tot_i);
    fftw_free(FFTout_Vx_tot_e);
    fftw_free(FFTout_Vy_tot_e);
    fftw_free(FFTout_Vz_tot_e);
    fftw_free(FFTout_Vx_tot_i);
    fftw_free(FFTout_Vy_tot_i);
    fftw_free(FFTout_Vz_tot_i);
    fftw_free(FFTout_B);
    fftw_free(FFTout_E);
    fftw_free(FFTout_V_tot_e);
    fftw_free(FFTout_V_tot_i);

    fftw_destroy_plan(plan_Bx);
    fftw_destroy_plan(plan_By);
    fftw_destroy_plan(plan_Bz);
    fftw_destroy_plan(plan_Ex);
    fftw_destroy_plan(plan_Ey);
    fftw_destroy_plan(plan_Ez);
    fftw_destroy_plan(plan_Jx_tot_e);
    fftw_destroy_plan(plan_Jy_tot_e);
    fftw_destroy_plan(plan_Jz_tot_e);
    fftw_destroy_plan(plan_Jx_tot_i);
    fftw_destroy_plan(plan_Jy_tot_i);
    fftw_destroy_plan(plan_Jz_tot_i);
    fftw_destroy_plan(plan_rho_tot_e);
    fftw_destroy_plan(plan_rho_tot_i);
    fftw_destroy_plan(plan_Vx_tot_e);
    fftw_destroy_plan(plan_Vy_tot_e);
    fftw_destroy_plan(plan_Vz_tot_e);
    fftw_destroy_plan(plan_Vx_tot_i);
    fftw_destroy_plan(plan_Vy_tot_i);
    fftw_destroy_plan(plan_Vz_tot_i);
    fftw_destroy_plan(plan_B);
    fftw_destroy_plan(plan_E);
    fftw_destroy_plan(plan_V_tot_e);
    fftw_destroy_plan(plan_V_tot_i);

    return 0;
}
