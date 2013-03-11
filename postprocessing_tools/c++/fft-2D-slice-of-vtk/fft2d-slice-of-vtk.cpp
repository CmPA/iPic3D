/***************************************************************************                        
    Extract 2D slice from a 3D VTK from iPIC3D and perform a 2D FFT in space

    Author               : Alexander E. Vapirev, CPA, KU Leuven

    Begin                : 18 Oct 2011

    Help                 : This c++ code extracts a user defined fragment
                           of any quantity from an iPIC3D style VTK and then perform and 2D FFT
                           in space of that slice fragment. The VTKs must be written
                           using the hdf-to-vtk-converter.cpp by A.E.Vapirev which is based on
                           the original HDFtoVTK.cpp by S.Markidis.

    Output               : The absolute value of the FFT is being written out and the x-axis is frequency
                           in radians/wci which are the frequencies correspoinding to the FFT elements.
                           The y-axis (spatial dimension), i.e., the line of virtual probes in space
                           is wave number in k/d_i which are the wave numbers corresponding to the FFT elements.
                           Then the frequency and wavenumbers are scaled - see the Note below.

    Change log           :

    If using FFTW3, compile with "c++ fft-virtsat-point.cpp -I/usr/local/include/ -L/usr/lib/ -lfftw3 -lm"
                         or with "c++ fft-virtsat-point.cpp -lfftw3 -lm"
                         of if you install FFTW3 from source adjust the 'include' and 'lib' locations accordingly.
                         See here for examples http://opencv-code.com/Using_the_Fast_Fourier_Transform_Library_FFTW
                         
    Note:                The FFT contains information between 0 and fs, where fs is the sample frequency.
                         The sampling frequency is the frequency rate between the adjacent samples.
                         It is calculated as fs=1/delta_T. In the case of iPIC3D delta_T=wci*dt
                         So if the frequency axis is plotted as (0:(M-1))/M it is actually f/fs what is plotted over X-axis.
                         If we multiply by 2Pi/delta_T, e.g. [(0:(M-1))/M]*[2Pi/delta_T] we get it in radians/wci.
                         If we multiply by 2Pi/dt, e.g. [(0:(M-1))/M]*[2Pi/dt] we get it in radians.
                         
                         For the y-axis the information is between 0 and ks, where ks is the sampling wave number.
                         The sampling wave number is the wave numbers corresponding to the wavenumber between adjacent
                         probes and is calculated as ks=1/delta_distance_between_probes. In the case of iPIC3D
                         delta_distance_between_probes has to be looked up for each separate run for each direction.
                         Most probably it will be simulation box size divided by the (number of probes - 1) in that dimension.
                         So if the discrete wave numbers are plotted as (0:(N-1))/N it is actually k/d_i what
                         is plotted over y-axis.
                         If we multuply by 2Pi/delta_Y, e.g., [(0:(N-1))/N]*[2Pi/delta_Y] we get it in k/d_i.
                         
                         All this scaling should be done in the plotting program. This script only computes values.


 ************************************************************************** */

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

using std::string;
using std::stringstream;
using std::ofstream;
using std::cout;
using std::endl;

int main(int argc, char* argv[]) {
/*
    // example input in file "user_inputs.h"
    char vtk_filename[]="rho_tot_electrons_cycle10000.vtk";
    string slice_normal="y";
    // var_type can be "scalar" or "vectorX", "vectorY", "vectorZ" for extracting a vector component
    string var_type="scalar";
    float slice_pos=7.5;
    // begin and end coordinates for the portion of the slice in 2D
    float coord11=0.0;
    float coord12=40.0;
    float coord21=0.0;
    float coord22=15.0;
*/

    #include "user_inputs.h"

    slice_pos = slice_pos - 1.0;    // c++ adds "1" so we need to subtract "1"

    int i,j,k,ii,jj,kk,idummy;
    char cdummy[]="";
    float fdummy;

    // index dimensions in x,y,z
    int nx,ny,nz;

    // origin of the 3D box
    float x0,y0,z0;

    // spacing in x,y,z
    float dx,dy,dz;

    // begin and end indices for the portion of the slice in 2D
    int index11;
    int index12;
    int index21;
    int index22;

    // open the VTK file
    FILE* vtk_file=fopen(vtk_filename,"r");

    // read the VTK header - dimensions, spacing, etc...
    cout << "---> In VTK file " << vtk_filename << ":" << endl;

    // skip (void read) first four lines of header in file
    for (k=0;k<4;k++){
        do idummy = fgetc(vtk_file); while (idummy != '\n');
    }

    // read dimensions
    fscanf(vtk_file,"%s %i %i %i",cdummy,&nx,&ny,&nz);
    cout << "---> Got dimensions: " << nx << " " << ny << " " << nz << endl;

    // read origin
    fscanf(vtk_file,"%s %f %f %f",cdummy,&x0,&y0,&z0);
    cout << "---> Got origin:     " << x0 << " " << y0 << " " << z0 << endl;

    // read spacing for the structured grid
    fscanf(vtk_file,"%s %f %f %f",cdummy,&dx,&dy,&dz);
    cout << "---> Got spacing:    " << dx << " " << dy << " " << dz << endl;

    // skip (void read) the rest of header in file
    for (k=0;k<5;k++){
        do idummy = fgetc(vtk_file); while (idummy != '\n');
    }

    // finished reading the VTK header

    cout << "---> Determining array indices corresponding to the required spatial subset in the 2D slice:" << endl;

    float closest_slice_pos;
    int index_closest_slice_pos;

    // find the index corresponding to the closest cut to slice_pos in x, y, or z

    if ( slice_normal == "x" ){
       closest_slice_pos=nx*dx;
       for (ii=0; ii < nx; ii++){
           if ( abs(slice_pos-ii*dx) <= closest_slice_pos ){
              closest_slice_pos=abs(slice_pos-ii*dx);
              index_closest_slice_pos=ii;
           }
       }
    }

    if ( slice_normal == "y" ){
       closest_slice_pos=ny*dy;
       for (jj=0; jj < ny; jj++){
           if ( abs(slice_pos-jj*dy) <= closest_slice_pos ){
              closest_slice_pos=abs(slice_pos-jj*dy);
              index_closest_slice_pos=jj;
           }
       }
    }

    if ( slice_normal == "z" ){
       closest_slice_pos=nz*dz;
       for (kk=0; kk < nz; kk++){
           if ( abs(slice_pos-kk*dz) <= closest_slice_pos ){
              closest_slice_pos=abs(slice_pos-kk*dz);
              index_closest_slice_pos=kk;
           }
       }
    }

    cout << "---> Desired 2D slice at " << slice_normal << "=" << slice_pos + 1.0 << endl;
    cout << "---> Found closest position difference " << closest_slice_pos << endl;
    cout << "---> Found index = " << index_closest_slice_pos << " corresponding to the closest slice position along the " << slice_normal << "-normal" << endl;

    // Find the coefficients corresponding to the desired dimensional beginning and ending of the slice

    if ( slice_normal == "x" ){
        for (jj=0; jj < ny;jj++) if (jj*dy <= coord11) index11=jj;
        for (jj=0; jj < ny;jj++) if (jj*dy <= coord12) index12=jj;
        for (kk=0; kk < nz;kk++) if (kk*dz <= coord21) index21=kk;
        for (kk=0; kk < nz;kk++) if (kk*dz <= coord22) index22=kk;
    }

    if ( slice_normal == "y" ){
        for (ii=0; ii < nx;ii++) if (ii*dx <= coord11) index11=ii;
        for (ii=0; ii < nx;ii++) if (ii*dx <= coord12) index12=ii;
        for (kk=0; kk < nz;kk++) if (kk*dz <= coord21) index21=kk;
        for (kk=0; kk < nz;kk++) if (kk*dz <= coord22) index22=kk;
    }

    if ( slice_normal == "z" ){
        for (ii=0; ii < nx;ii++) if (ii*dx <= coord11) index11=ii;
        for (ii=0; ii < nx;ii++) if (ii*dx <= coord12) index12=ii;
        for (jj=0; jj < ny;jj++) if (jj*dy <= coord21) index21=jj;
        for (jj=0; jj < ny;jj++) if (jj*dy <= coord22) index22=jj;
    }

    cout << "---> Found begin and end indices for dimension 1 and 2:" << endl;
    cout << "---> index11 = " << index11 << " | index12 = " << index12 << endl;
    cout << "---> index21 = " << index21 << " | index22 = " << index22 << endl;

    // write the indices in a file
    FILE* file_indices=fopen("indices.txt","w");
    fprintf(file_indices,"%i %i %i %i",index11,index12,index21,index22);
    fclose(file_indices);

    // define the subset which we will extract
    const int n1 = index12 - index11;
    const int n2 = index22 - index21;
    float scalar[n1+1][n2+1];

    // read the data

    cout << "---> Extracting the slice fragment from VTK..." << endl;

    float var_value,var_value_X,var_value_Y,var_value_Z;

    if ( slice_normal == "x" ){
       for (kk=0; kk < nz;kk++){
           for (jj=0; jj < ny;jj++){
               for (ii=0; ii < nx;ii++){
                   if (var_type == "scalar" ) fscanf(vtk_file,"%f",&var_value);                        // reading scalar
                   if (var_type == "vectorX" ) fscanf(vtk_file,"%f %f %f",&var_value,&fdummy,&fdummy); // reading x-component
                   if (var_type == "vectorY" ) fscanf(vtk_file,"%f %f %f",&fdummy,&var_value,&fdummy); // reading y-component
                   if (var_type == "vectorZ" ) fscanf(vtk_file,"%f %f %f",&fdummy,&fdummy,&var_value); // reading z-component
                   if (var_type == "vector_magnitude" ){                                               // reading vector magnitude
                      fscanf(vtk_file,"%f %f %f",&var_value_X,&var_value_Y,&var_value_Z);
                      var_value=sqrt(var_value_X*var_value_X+var_value_Y*var_value_Y+var_value_Z*var_value_Z);
                   }
                   // read only the fraction limited by indices (11,12),(21,22) at index_closest_slice_pos
                   if (ii == index_closest_slice_pos && jj >= index11 && jj <= index12 && kk >= index21 && kk <= index22){
                      scalar[jj-index11][kk-index21]=var_value;
                   }
               }
           }
       }
    }
    
    if ( slice_normal == "y" ){
       for (kk=0; kk < nz;kk++){
           for (jj=0; jj < ny;jj++){
               for (ii=0; ii < nx;ii++){
                   if (var_type == "scalar" ) fscanf(vtk_file,"%f",&var_value);                        // reading scalar
                   if (var_type == "vectorX" ) fscanf(vtk_file,"%f %f %f",&var_value,&fdummy,&fdummy); // reading x-component
                   if (var_type == "vectorY" ) fscanf(vtk_file,"%f %f %f",&fdummy,&var_value,&fdummy); // reading y-component
                   if (var_type == "vectorZ" ) fscanf(vtk_file,"%f %f %f",&fdummy,&fdummy,&var_value); // reading z-component
                   if (var_type == "vector_magnitude" ){                                               // reading vector magnitude
                      fscanf(vtk_file,"%f %f %f",&var_value_X,&var_value_Y,&var_value_Z);
                      var_value=sqrt(var_value_X*var_value_X+var_value_Y*var_value_Y+var_value_Z*var_value_Z);
                   }
                   // read only the fraction limited by indices (11,12),(21,22) at index_closest_slice_pos
                   if (jj == index_closest_slice_pos && ii >= index11 && ii <= index12 && kk >= index21 && kk <= index22){
                      scalar[ii-index11][kk-index21]=var_value;
                   }
               }
           }
       }
    }

    if ( slice_normal == "z" ){
       for (kk=0; kk < nz;kk++){
           for (jj=0; jj < ny;jj++){
               for (ii=0; ii < nx;ii++){
                   if (var_type == "scalar" ) fscanf(vtk_file,"%f",&var_value);                        // reading scalar
                   if (var_type == "vectorX" ) fscanf(vtk_file,"%f %f %f",&var_value,&fdummy,&fdummy); // reading x-component
                   if (var_type == "vectorY" ) fscanf(vtk_file,"%f %f %f",&fdummy,&var_value,&fdummy); // reading y-component
                   if (var_type == "vectorZ" ) fscanf(vtk_file,"%f %f %f",&fdummy,&fdummy,&var_value); // reading z-component
                   if (var_type == "vector_magnitude" ){                                               // reading vector magnitude
                      fscanf(vtk_file,"%f %f %f",&var_value_X,&var_value_Y,&var_value_Z);
                      var_value=sqrt(var_value_X*var_value_X+var_value_Y*var_value_Y+var_value_Z*var_value_Z);
                   }
                   // read only the fraction limited by indices (11,12),(21,22) at index_closest_slice_pos
                   if (kk == index_closest_slice_pos && ii >= index11 && ii <= index12 && jj >= index21 && jj <= index22){
                      scalar[ii-index11][jj-index21]=var_value;
                   }
               }
           }
       }
    }
    
    fclose(vtk_file);

    // perform the 2D FFT

    // convert 2D array to 1D for FFT conversion

    int NN = (n1+1)*(n2+1);
    double scalar_2DFFT[NN];
    int inn = 0;

    for (i=0; i<n1+1; i++){
        for (j=0; j<n2+1; j++){
            scalar_2DFFT[inn]=scalar[i][j];
            inn = inn + 1;
        }
    }

    // Allocate array for DFT result

    int NNfft2d = (n1+1);
    int MMfft2d = (n2+1);
    int MM = NNfft2d*(MMfft2d/2+1);
    fftw_complex* forward_2DFFT;
    forward_2DFFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * MM);

    // define plan

    fftw_plan plan_forward_2DFFT;
    plan_forward_2DFFT = fftw_plan_dft_r2c_2d(NNfft2d, MMfft2d, scalar_2DFFT, forward_2DFFT, FFTW_ESTIMATE);

    // Compute forward discrete FT in 2D
    fftw_execute(plan_forward_2DFFT);

    // write out for gnuplot - the original quantity and its FFT - column number is x and row number is y

    cout << "---> Writing out the slice fragment and its 2D FFT..." << endl;

    // write out the extracted slice fragment

    FILE* output_file=fopen("output_2D_slice_VTK.txt","w");

    for (j=0; j<MMfft2d; j++){
        for (i=0; i<NNfft2d; i++){
            fprintf(output_file,"%E ",scalar[i][j]);
        }
        fprintf(output_file,"\n","");
    }    

    fclose(output_file);

    // output the 2D FFT of the slice fragment

    FILE* output_file_FFT=fopen("output_fft_2D_slice_VTK.txt","w");
    float absFFT2D;
    float forward_2DFFT_array2D_real[NNfft2d][MMfft2d/2+1];
    float forward_2DFFT_array2D_imaginary[NNfft2d][MMfft2d/2+1];

    // first convert the row-major output into a 2D array

    int ifft=0;
    for (i=0; i<NNfft2d; i++){
        for (j=0; j<(MMfft2d/2+1); j++){
            forward_2DFFT_array2D_real[i][j]      = forward_2DFFT[ifft][0];
            forward_2DFFT_array2D_imaginary[i][j] = forward_2DFFT[ifft][1];
            ifft=ifft+1;
        }
    }    

    //  and then write it out

    for (j=0; j<(MMfft2d/2+1); j++){
        for (i=0; i<NNfft2d; i++){
            absFFT2D=sqrt( forward_2DFFT_array2D_real[i][j]*forward_2DFFT_array2D_real[i][j] + forward_2DFFT_array2D_imaginary[i][j]*forward_2DFFT_array2D_imaginary[i][j] );
            fprintf(output_file_FFT,"%E ",absFFT2D);
        }
        fprintf(output_file_FFT,"\n","");
    }    

    fclose(output_file_FFT);

    // empty the memory
    fftw_destroy_plan(plan_forward_2DFFT);
    fftw_free(forward_2DFFT);

    return(0);
}   // end main program

