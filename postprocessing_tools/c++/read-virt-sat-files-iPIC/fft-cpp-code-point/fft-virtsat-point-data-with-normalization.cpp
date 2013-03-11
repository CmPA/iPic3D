/***************************************************************************

    fft-virtsat-point-data.cpp  -  FFT of 1D real data extracted with process-virtual-sat-point.cpp

    Author:  A.E.Vapirev, KU Leuven, Afdeling Plasma-astrofysica
             2011 Oct 7

    Output               : The absolute value of the FFT is being written out and the x-axis is labeled
                           in radians as w = 2*pi*(0:(N-1))/N which are the frequencies correspoinding to the FFT elements.
                           (The length of the descrete FT corresponds to the frequency 2*Pi in the descrete time FT.)

    Change log           :

    If using FFTW3, compile with "c++ fft-virtsat-point.cpp -I/usr/local/include/ -L/usr/lib/ -lfftw3 -lm"
                         or with "c++ fft-virtsat-point.cpp -lfftw3 -lm"
                         of if you install FFTW3 from source adjust the 'include' and 'lib' locations accordingly.
                         See here for examples http://opencv-code.com/Using_the_Fast_Fourier_Transform_Library_FFTW

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

//    sscanf(argv[1],"%d",&number_of_data_rows);

    // get the number of data points from 'BestVirtualSatelliteInfoFoundPoint.txt'
    char VirtualSatelliteFilenameFound[] = "";  // dummy
    int satellite_number_best_match;            // dummy
    int number_of_data_rows;
    char info_fname[] = "BestVirtualSatelliteInfoFoundPoint.txt";
    FILE *satellite_info_file;
    satellite_info_file = fopen(info_fname,"r");
    fscanf(satellite_info_file,"%s %i %i",VirtualSatelliteFilenameFound,&satellite_number_best_match,&number_of_data_rows);
    fclose(satellite_info_file);

    int N=number_of_data_rows;   // number of input data points
    int M=N/2+1;                 // number of output FFT data points

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
    double * FFTin_Bx  = (double*) malloc(sizeof(double) * N);
    double * FFTin_By  = (double*) malloc(sizeof(double) * N);
    double * FFTin_Bz  = (double*) malloc(sizeof(double) * N);
    double * FFTin_Ex  = (double*) malloc(sizeof(double) * N);
    double * FFTin_Ey  = (double*) malloc(sizeof(double) * N);
    double * FFTin_Ez  = (double*) malloc(sizeof(double) * N);
    double * FFTin_Jx_tot_e  = (double*) malloc(sizeof(double) * N);
    double * FFTin_Jy_tot_e  = (double*) malloc(sizeof(double) * N);
    double * FFTin_Jz_tot_e  = (double*) malloc(sizeof(double) * N);
    double * FFTin_Jx_tot_i  = (double*) malloc(sizeof(double) * N);
    double * FFTin_Jy_tot_i  = (double*) malloc(sizeof(double) * N);
    double * FFTin_Jz_tot_i  = (double*) malloc(sizeof(double) * N);
    double * FFTin_rho_tot_e  = (double*) malloc(sizeof(double) * N);
    double * FFTin_rho_tot_i  = (double*) malloc(sizeof(double) * N);
    double * FFTin_Vx_tot_e  = (double*) malloc(sizeof(double) * N);
    double * FFTin_Vy_tot_e  = (double*) malloc(sizeof(double) * N);
    double * FFTin_Vz_tot_e  = (double*) malloc(sizeof(double) * N);
    double * FFTin_Vx_tot_i  = (double*) malloc(sizeof(double) * N);
    double * FFTin_Vy_tot_i  = (double*) malloc(sizeof(double) * N);
    double * FFTin_Vz_tot_i  = (double*) malloc(sizeof(double) * N);
    double * FFTin_B  = (double*) malloc(sizeof(double) * N);
    double * FFTin_E  = (double*) malloc(sizeof(double) * N);
    double * FFTin_V_tot_e  = (double*) malloc(sizeof(double) * N);
    double * FFTin_V_tot_i  = (double*) malloc(sizeof(double) * N);
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
    fftw_plan plan_Bx = fftw_plan_dft_r2c_1d(N, FFTin_Bx, FFTout_Bx, FFTW_ESTIMATE);
    fftw_plan plan_By = fftw_plan_dft_r2c_1d(N, FFTin_By, FFTout_By, FFTW_ESTIMATE);
    fftw_plan plan_Bz = fftw_plan_dft_r2c_1d(N, FFTin_Bz, FFTout_Bz, FFTW_ESTIMATE);
    fftw_plan plan_Ex = fftw_plan_dft_r2c_1d(N, FFTin_Ex, FFTout_Ex, FFTW_ESTIMATE);
    fftw_plan plan_Ey = fftw_plan_dft_r2c_1d(N, FFTin_Ey, FFTout_Ey, FFTW_ESTIMATE);
    fftw_plan plan_Ez = fftw_plan_dft_r2c_1d(N, FFTin_Ez, FFTout_Ez, FFTW_ESTIMATE);
    fftw_plan plan_Jx_tot_e = fftw_plan_dft_r2c_1d(N, FFTin_Jx_tot_e, FFTout_Jx_tot_e, FFTW_ESTIMATE);
    fftw_plan plan_Jy_tot_e = fftw_plan_dft_r2c_1d(N, FFTin_Jy_tot_e, FFTout_Jy_tot_e, FFTW_ESTIMATE);
    fftw_plan plan_Jz_tot_e = fftw_plan_dft_r2c_1d(N, FFTin_Jz_tot_e, FFTout_Jz_tot_e, FFTW_ESTIMATE);
    fftw_plan plan_Jx_tot_i = fftw_plan_dft_r2c_1d(N, FFTin_Jx_tot_i, FFTout_Jx_tot_i, FFTW_ESTIMATE);
    fftw_plan plan_Jy_tot_i = fftw_plan_dft_r2c_1d(N, FFTin_Jy_tot_i, FFTout_Jy_tot_i, FFTW_ESTIMATE);
    fftw_plan plan_Jz_tot_i = fftw_plan_dft_r2c_1d(N, FFTin_Jz_tot_i, FFTout_Jz_tot_i, FFTW_ESTIMATE);
    fftw_plan plan_rho_tot_e = fftw_plan_dft_r2c_1d(N, FFTin_rho_tot_e, FFTout_rho_tot_e, FFTW_ESTIMATE);
    fftw_plan plan_rho_tot_i = fftw_plan_dft_r2c_1d(N, FFTin_rho_tot_i, FFTout_rho_tot_i, FFTW_ESTIMATE);
    fftw_plan plan_Vx_tot_e = fftw_plan_dft_r2c_1d(N, FFTin_Vx_tot_e, FFTout_Vx_tot_e, FFTW_ESTIMATE);
    fftw_plan plan_Vy_tot_e = fftw_plan_dft_r2c_1d(N, FFTin_Vy_tot_e, FFTout_Vy_tot_e, FFTW_ESTIMATE);
    fftw_plan plan_Vz_tot_e = fftw_plan_dft_r2c_1d(N, FFTin_Vz_tot_e, FFTout_Vz_tot_e, FFTW_ESTIMATE);
    fftw_plan plan_Vx_tot_i = fftw_plan_dft_r2c_1d(N, FFTin_Vx_tot_i, FFTout_Vx_tot_i, FFTW_ESTIMATE);
    fftw_plan plan_Vy_tot_i = fftw_plan_dft_r2c_1d(N, FFTin_Vy_tot_i, FFTout_Vy_tot_i, FFTW_ESTIMATE);
    fftw_plan plan_Vz_tot_i = fftw_plan_dft_r2c_1d(N, FFTin_Vz_tot_i, FFTout_Vz_tot_i, FFTW_ESTIMATE);
    fftw_plan plan_B = fftw_plan_dft_r2c_1d(N, FFTin_B, FFTout_B, FFTW_ESTIMATE);
    fftw_plan plan_E = fftw_plan_dft_r2c_1d(N, FFTin_E, FFTout_E, FFTW_ESTIMATE);
    fftw_plan plan_V_tot_e = fftw_plan_dft_r2c_1d(N, FFTin_V_tot_e, FFTout_V_tot_e, FFTW_ESTIMATE);
    fftw_plan plan_V_tot_i = fftw_plan_dft_r2c_1d(N, FFTin_V_tot_i, FFTout_V_tot_i, FFTW_ESTIMATE);

    // read from input data file
    char data_fname[] = "sat_output_point.txt";
    FILE *input_file;
    input_file = fopen(data_fname,"r");

    // skip (void read) first line (header) in file
    int idummy;
    do idummy = fgetc(input_file); while (idummy != '\n');

    for (int i=0;i<number_of_data_rows;i++){                                // loop time records
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
        FFTin_Bx[i] = Bx[i];
        FFTin_By[i] = By[i];
        FFTin_Bz[i] = Bz[i];
        FFTin_Ex[i] = Ex[i];
        FFTin_Ey[i] = Ey[i];
        FFTin_Ez[i] = Ez[i];
        FFTin_Jx_tot_e[i] = Jx_tot_e[i];
        FFTin_Jy_tot_e[i] = Jy_tot_e[i];
        FFTin_Jz_tot_e[i] = Jz_tot_e[i];
        FFTin_Jx_tot_i[i] = Jx_tot_i[i];
        FFTin_Jy_tot_i[i] = Jy_tot_i[i];
        FFTin_Jz_tot_i[i] = Jz_tot_i[i];
        FFTin_rho_tot_e[i] = rho_tot_e[i];
        FFTin_rho_tot_i[i] = rho_tot_i[i];
        FFTin_Vx_tot_e[i] = Vx_tot_e[i];
        FFTin_Vy_tot_e[i] = Vy_tot_e[i];
        FFTin_Vz_tot_e[i] = Vz_tot_e[i];
        FFTin_Vx_tot_i[i] = Vx_tot_i[i];
        FFTin_Vy_tot_i[i] = Vy_tot_i[i];
        FFTin_Vz_tot_i[i] = Vz_tot_i[i];
        FFTin_B[i] = B[i];
        FFTin_E[i] = E[i];
        FFTin_V_tot_e[i] = V_tot_e[i];
        FFTin_V_tot_i[i] = V_tot_i[i];
    }

    // close input file
    fclose(input_file);

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

    // define max of abs values of the FFT output
    float max_absFFTout_Bx=1.0e-30;
    float max_absFFTout_By=1.0e-30;
    float max_absFFTout_Bz=1.0e-30;
    float max_absFFTout_Ex=1.0e-30;
    float max_absFFTout_Ey=1.0e-30;
    float max_absFFTout_Ez=1.0e-30;
    float max_absFFTout_Jx_tot_e=1.0e-30;
    float max_absFFTout_Jy_tot_e=1.0e-30;
    float max_absFFTout_Jz_tot_e=1.0e-30;
    float max_absFFTout_Jx_tot_i=1.0e-30;
    float max_absFFTout_Jy_tot_i=1.0e-30;
    float max_absFFTout_Jz_tot_i=1.0e-30;
    float max_absFFTout_rho_tot_e=1.0e-30;
    float max_absFFTout_rho_tot_i=1.0e-30;
    float max_absFFTout_Vx_tot_e=1.0e-30;
    float max_absFFTout_Vy_tot_e=1.0e-30;
    float max_absFFTout_Vz_tot_e=1.0e-30;
    float max_absFFTout_Vx_tot_i=1.0e-30;
    float max_absFFTout_Vy_tot_i=1.0e-30;
    float max_absFFTout_Vz_tot_i=1.0e-30;
    float max_absFFTout_B=1.0e-30;
    float max_absFFTout_E=1.0e-30;
    float max_absFFTout_V_tot_e=1.0e-30;
    float max_absFFTout_V_tot_i=1.0e-30;

    // find max of abs values of the FFT output
    for (int i=0;i<M;i++){                                // loop time records
        absFFTout_Bx=sqrt(FFTout_Bx[i][0]*FFTout_Bx[i][0]+FFTout_Bx[i][1]*FFTout_Bx[i][1]);
        if (absFFTout_Bx>max_absFFTout_Bx) max_absFFTout_Bx=absFFTout_Bx;
        absFFTout_By=sqrt(FFTout_By[i][0]*FFTout_By[i][0]+FFTout_By[i][1]*FFTout_By[i][1]);
        if (absFFTout_By>max_absFFTout_By) max_absFFTout_By=absFFTout_By;
        absFFTout_Bz=sqrt(FFTout_Bz[i][0]*FFTout_Bz[i][0]+FFTout_Bz[i][1]*FFTout_Bz[i][1]);
        if (absFFTout_Bz>max_absFFTout_Bz) max_absFFTout_Bz=absFFTout_Bz;

        absFFTout_Ex=sqrt(FFTout_Ex[i][0]*FFTout_Ex[i][0]+FFTout_Ex[i][1]*FFTout_Ex[i][1]);
        if (absFFTout_Ex>max_absFFTout_Ex) max_absFFTout_Ex=absFFTout_Ex;
        absFFTout_Ey=sqrt(FFTout_Ey[i][0]*FFTout_Ey[i][0]+FFTout_Ey[i][1]*FFTout_Ey[i][1]);
        if (absFFTout_Ey>max_absFFTout_Ey) max_absFFTout_Ey=absFFTout_Ey;
        absFFTout_Ez=sqrt(FFTout_Ez[i][0]*FFTout_Ez[i][0]+FFTout_Ez[i][1]*FFTout_Ez[i][1]);
        if (absFFTout_Ez>max_absFFTout_Ez) max_absFFTout_Ez=absFFTout_Ez;

        absFFTout_Jx_tot_e=sqrt(FFTout_Jx_tot_e[i][0]*FFTout_Jx_tot_e[i][0]+FFTout_Jx_tot_e[i][1]*FFTout_Jx_tot_e[i][1]);
        if (absFFTout_Jx_tot_e>max_absFFTout_Jx_tot_e) max_absFFTout_Jx_tot_e=absFFTout_Jx_tot_e;
        absFFTout_Jy_tot_e=sqrt(FFTout_Jy_tot_e[i][0]*FFTout_Jy_tot_e[i][0]+FFTout_Jy_tot_e[i][1]*FFTout_Jy_tot_e[i][1]);
        if (absFFTout_Jy_tot_e>max_absFFTout_Jy_tot_e) max_absFFTout_Jy_tot_e=absFFTout_Jy_tot_e;
        absFFTout_Jz_tot_e=sqrt(FFTout_Jz_tot_e[i][0]*FFTout_Jz_tot_e[i][0]+FFTout_Jz_tot_e[i][1]*FFTout_Jz_tot_e[i][1]);
        if (absFFTout_Jz_tot_e>max_absFFTout_Jz_tot_e) max_absFFTout_Jz_tot_e=absFFTout_Jz_tot_e;

        absFFTout_Jx_tot_i=sqrt(FFTout_Jx_tot_i[i][0]*FFTout_Jx_tot_i[i][0]+FFTout_Jx_tot_i[i][1]*FFTout_Jx_tot_i[i][1]);
        if (absFFTout_Jx_tot_i>max_absFFTout_Jx_tot_i) max_absFFTout_Jx_tot_i=absFFTout_Jx_tot_i;
        absFFTout_Jy_tot_i=sqrt(FFTout_Jy_tot_i[i][0]*FFTout_Jy_tot_i[i][0]+FFTout_Jy_tot_i[i][1]*FFTout_Jy_tot_i[i][1]);
        if (absFFTout_Jy_tot_i>max_absFFTout_Jy_tot_i) max_absFFTout_Jy_tot_i=absFFTout_Jy_tot_i;
        absFFTout_Jz_tot_i=sqrt(FFTout_Jz_tot_i[i][0]*FFTout_Jz_tot_i[i][0]+FFTout_Jz_tot_i[i][1]*FFTout_Jz_tot_i[i][1]);
        if (absFFTout_Jz_tot_i>max_absFFTout_Jz_tot_i) max_absFFTout_Jz_tot_i=absFFTout_Jz_tot_i;

        absFFTout_rho_tot_e=sqrt(FFTout_rho_tot_e[i][0]*FFTout_rho_tot_e[i][0]+FFTout_rho_tot_e[i][1]*FFTout_rho_tot_e[i][1]);
        if (absFFTout_rho_tot_e>max_absFFTout_rho_tot_e) max_absFFTout_rho_tot_e=absFFTout_rho_tot_e;
        absFFTout_rho_tot_i=sqrt(FFTout_rho_tot_i[i][0]*FFTout_rho_tot_i[i][0]+FFTout_rho_tot_i[i][1]*FFTout_rho_tot_i[i][1]);
        if (absFFTout_rho_tot_i>max_absFFTout_rho_tot_i) max_absFFTout_rho_tot_i=absFFTout_rho_tot_i;

        absFFTout_Vx_tot_e=sqrt(FFTout_Vx_tot_e[i][0]*FFTout_Vx_tot_e[i][0]+FFTout_Vx_tot_e[i][1]*FFTout_Vx_tot_e[i][1]);
        if (absFFTout_Vx_tot_e>max_absFFTout_Vx_tot_e) max_absFFTout_Vx_tot_e=absFFTout_Vx_tot_e;
        absFFTout_Vy_tot_e=sqrt(FFTout_Vy_tot_e[i][0]*FFTout_Vy_tot_e[i][0]+FFTout_Vy_tot_e[i][1]*FFTout_Vy_tot_e[i][1]);
        if (absFFTout_Vy_tot_e>max_absFFTout_Vy_tot_e) max_absFFTout_Vy_tot_e=absFFTout_Vy_tot_e;
        absFFTout_Vz_tot_e=sqrt(FFTout_Vz_tot_e[i][0]*FFTout_Vz_tot_e[i][0]+FFTout_Vz_tot_e[i][1]*FFTout_Vz_tot_e[i][1]);
        if (absFFTout_Vz_tot_e>max_absFFTout_Vz_tot_e) max_absFFTout_Vz_tot_e=absFFTout_Vz_tot_e;

        absFFTout_Vx_tot_i=sqrt(FFTout_Vx_tot_i[i][0]*FFTout_Vx_tot_i[i][0]+FFTout_Vx_tot_i[i][1]*FFTout_Vx_tot_i[i][1]);
        if (absFFTout_Vx_tot_i>max_absFFTout_Vx_tot_i) max_absFFTout_Vx_tot_i=absFFTout_Vx_tot_i;
        absFFTout_Vy_tot_i=sqrt(FFTout_Vy_tot_i[i][0]*FFTout_Vy_tot_i[i][0]+FFTout_Vy_tot_i[i][1]*FFTout_Vy_tot_i[i][1]);
        if (absFFTout_Vy_tot_i>max_absFFTout_Vy_tot_i) max_absFFTout_Vy_tot_i=absFFTout_Vy_tot_i;
        absFFTout_Vz_tot_i=sqrt(FFTout_Vz_tot_i[i][0]*FFTout_Vz_tot_i[i][0]+FFTout_Vz_tot_i[i][1]*FFTout_Vz_tot_i[i][1]);
        if (absFFTout_Vz_tot_i>max_absFFTout_Vz_tot_i) max_absFFTout_Vz_tot_i=absFFTout_Vz_tot_i;

        absFFTout_B=sqrt(FFTout_B[i][0]*FFTout_B[i][0]+FFTout_B[i][1]*FFTout_B[i][1]);
        if (absFFTout_B>max_absFFTout_B) max_absFFTout_B=absFFTout_B;
        absFFTout_E=sqrt(FFTout_E[i][0]*FFTout_E[i][0]+FFTout_E[i][1]*FFTout_E[i][1]);
        if (absFFTout_E>max_absFFTout_E) max_absFFTout_E=absFFTout_E;

        absFFTout_V_tot_e=sqrt(FFTout_V_tot_e[i][0]*FFTout_V_tot_e[i][0]+FFTout_V_tot_e[i][1]*FFTout_V_tot_e[i][1]);
        if (absFFTout_V_tot_e>max_absFFTout_V_tot_e) max_absFFTout_V_tot_e=absFFTout_V_tot_e;
        absFFTout_V_tot_i=sqrt(FFTout_V_tot_i[i][0]*FFTout_V_tot_i[i][0]+FFTout_V_tot_i[i][1]*FFTout_V_tot_i[i][1]);
        if (absFFTout_V_tot_i>max_absFFTout_V_tot_i) max_absFFTout_V_tot_i=absFFTout_V_tot_i;
    }

    // the x-axis is frequency in radians
    float frequency_radians;

    // write out the data
    for (int i=0;i<M;i++){                                // loop time records
        frequency_radians=2.0*3.1415927*i/M;

        // find abs of FFT and normalize it
        absFFTout_Bx=sqrt(FFTout_Bx[i][0]*FFTout_Bx[i][0]+FFTout_Bx[i][1]*FFTout_Bx[i][1])/max_absFFTout_Bx;
        absFFTout_By=sqrt(FFTout_By[i][0]*FFTout_By[i][0]+FFTout_By[i][1]*FFTout_By[i][1])/max_absFFTout_By;
        absFFTout_Bz=sqrt(FFTout_Bz[i][0]*FFTout_Bz[i][0]+FFTout_Bz[i][1]*FFTout_Bz[i][1])/max_absFFTout_Bz;
        absFFTout_Ex=sqrt(FFTout_Ex[i][0]*FFTout_Ex[i][0]+FFTout_Ex[i][1]*FFTout_Ex[i][1])/max_absFFTout_Ez;
        absFFTout_Ey=sqrt(FFTout_Ey[i][0]*FFTout_Ey[i][0]+FFTout_Ey[i][1]*FFTout_Ey[i][1])/max_absFFTout_Ey;
        absFFTout_Ez=sqrt(FFTout_Ez[i][0]*FFTout_Ez[i][0]+FFTout_Ez[i][1]*FFTout_Ez[i][1])/max_absFFTout_Ez;
        absFFTout_Jx_tot_e=sqrt(FFTout_Jx_tot_e[i][0]*FFTout_Jx_tot_e[i][0]+FFTout_Jx_tot_e[i][1]*FFTout_Jx_tot_e[i][1])/max_absFFTout_Jx_tot_e;
        absFFTout_Jy_tot_e=sqrt(FFTout_Jy_tot_e[i][0]*FFTout_Jy_tot_e[i][0]+FFTout_Jy_tot_e[i][1]*FFTout_Jy_tot_e[i][1])/max_absFFTout_Jy_tot_e;
        absFFTout_Jz_tot_e=sqrt(FFTout_Jz_tot_e[i][0]*FFTout_Jz_tot_e[i][0]+FFTout_Jz_tot_e[i][1]*FFTout_Jz_tot_e[i][1])/max_absFFTout_Jz_tot_e;
        absFFTout_Jx_tot_i=sqrt(FFTout_Jx_tot_i[i][0]*FFTout_Jx_tot_i[i][0]+FFTout_Jx_tot_i[i][1]*FFTout_Jx_tot_i[i][1])/max_absFFTout_Jx_tot_i;
        absFFTout_Jy_tot_i=sqrt(FFTout_Jy_tot_i[i][0]*FFTout_Jy_tot_i[i][0]+FFTout_Jy_tot_i[i][1]*FFTout_Jy_tot_i[i][1])/max_absFFTout_Jy_tot_i;
        absFFTout_Jz_tot_i=sqrt(FFTout_Jz_tot_i[i][0]*FFTout_Jz_tot_i[i][0]+FFTout_Jz_tot_i[i][1]*FFTout_Jz_tot_i[i][1])/max_absFFTout_Jz_tot_i;
        absFFTout_rho_tot_e=sqrt(FFTout_rho_tot_e[i][0]*FFTout_rho_tot_e[i][0]+FFTout_rho_tot_e[i][1]*FFTout_rho_tot_e[i][1])/max_absFFTout_rho_tot_e;
        absFFTout_rho_tot_i=sqrt(FFTout_rho_tot_i[i][0]*FFTout_rho_tot_i[i][0]+FFTout_rho_tot_i[i][1]*FFTout_rho_tot_i[i][1])/max_absFFTout_rho_tot_i;
        absFFTout_Vx_tot_e=sqrt(FFTout_Vx_tot_e[i][0]*FFTout_Vx_tot_e[i][0]+FFTout_Vx_tot_e[i][1]*FFTout_Vx_tot_e[i][1])/max_absFFTout_Vx_tot_e;
        absFFTout_Vy_tot_e=sqrt(FFTout_Vy_tot_e[i][0]*FFTout_Vy_tot_e[i][0]+FFTout_Vy_tot_e[i][1]*FFTout_Vy_tot_e[i][1])/max_absFFTout_Vy_tot_e;
        absFFTout_Vz_tot_e=sqrt(FFTout_Vz_tot_e[i][0]*FFTout_Vz_tot_e[i][0]+FFTout_Vz_tot_e[i][1]*FFTout_Vz_tot_e[i][1])/max_absFFTout_Vz_tot_e;
        absFFTout_Vx_tot_i=sqrt(FFTout_Vx_tot_i[i][0]*FFTout_Vx_tot_i[i][0]+FFTout_Vx_tot_i[i][1]*FFTout_Vx_tot_i[i][1])/max_absFFTout_Vx_tot_i;
        absFFTout_Vy_tot_i=sqrt(FFTout_Vy_tot_i[i][0]*FFTout_Vy_tot_i[i][0]+FFTout_Vy_tot_i[i][1]*FFTout_Vy_tot_i[i][1])/max_absFFTout_Vy_tot_i;
        absFFTout_Vz_tot_i=sqrt(FFTout_Vz_tot_i[i][0]*FFTout_Vz_tot_i[i][0]+FFTout_Vz_tot_i[i][1]*FFTout_Vz_tot_i[i][1])/max_absFFTout_Vz_tot_i;
        absFFTout_B=sqrt(FFTout_B[i][0]*FFTout_B[i][0]+FFTout_B[i][1]*FFTout_B[i][1])/max_absFFTout_B;
        absFFTout_E=sqrt(FFTout_E[i][0]*FFTout_E[i][0]+FFTout_E[i][1]*FFTout_E[i][1])/max_absFFTout_E;
        absFFTout_V_tot_e=sqrt(FFTout_V_tot_e[i][0]*FFTout_V_tot_e[i][0]+FFTout_V_tot_e[i][1]*FFTout_V_tot_e[i][1])/max_absFFTout_V_tot_e;
        absFFTout_V_tot_i=sqrt(FFTout_V_tot_i[i][0]*FFTout_V_tot_i[i][0]+FFTout_V_tot_i[i][1]*FFTout_V_tot_i[i][1])/max_absFFTout_V_tot_i;

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
