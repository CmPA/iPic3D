/***************************************************************************

    process-virtual-sat.cpp  -  postprocessing program to extract data from VirtualSatellites* files from iPIC3D

    Author:  A.E.Vapirev, KU Leuven, Afdeling Plasma-astrofysica
             2012 May 22

    Change log           : this is a quick fix for a 2D version based on the 3D one and it has not been optimized at all

    If using FFTW3, compile with "c++ process-virtual-sat-point.cpp -I/usr/local/include/ -L/usr/lib/ -lfftw3 -lm"
                         or with "c++ process-virtual-sat-point.cpp -lfftw3 -lm"
                         of if you install FFTW3 from source adjust the 'include' and 'lib' locations accordingly.

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

using namespace std;

// FUNCTION DEFINITION

////////////////////////////////////////////////////////////////////////
// Function 'find_virtual_sat' finds the best satellite match with closest position to user given (x,y,z) spatial coordinates.
// The result is used by funciton 'read_virtual_sat'.
////////////////////////////////////////////////////////////////////////

void find_virtual_sat(float x_satellite_user,
                      float y_satellite_user,
                      int const nproc){

    FILE *virt_sat_file;                       // input file

    const int number_of_coords = 2;            // dimension
    const int number_of_satellites = 16;       // satellites per file

    int i,j,k,iproc;
    float temp_number;

    int satellite_number_best_match;                                        // later used in read-virtual-sat.cpp to extract the corresponding data
    string VirtualSatelliteFilenameFoundString;                             // this file contains the best satellite coordinates match

    vector<float> desired_sat_coordinates (number_of_coords,0);             // define 3D vector with zeros - the user desired satellite position
    vector<float> current_sat_coordinates (number_of_coords,0);             // define 3D vector with zeros - reading current coords. from the file
    vector<float> closest_sat_coordinates (number_of_coords,0);             // define 3D vector with zeros - closest satellite position to what we want
    vector<float> current_vector_difference (number_of_coords,0); 
    vector<float> closest_vector_difference (number_of_coords,0);           // this will end up being closest to the desired satellite position

    float mod_current_vector_difference;
    float mod_closest_vector_difference;

    // set up vector with the user defined coordinates
    desired_sat_coordinates[0] = x_satellite_user;
    desired_sat_coordinates[1] = y_satellite_user;

    // set up initially the closest satelite position to something much bigger than the desired position
    closest_vector_difference[0] = desired_sat_coordinates[0] * 1.0E30;
    closest_vector_difference[1] = desired_sat_coordinates[1] * 1.0E30;

    string temp_string;

    cout << "---> Searching for best satellite match..." << endl;

    // start best satellite search over all files
    for (int iproc=0; iproc < nproc; iproc++){                              // loop over all files
        stringstream ss;
        ss << iproc;
        temp_string = "VirtualSatelliteTraces" + ss.str() + ".txt";
        virt_sat_file = fopen(temp_string.c_str(),"r");
        if (virt_sat_file){                                                 // if the file exists
           for (i=0;i<number_of_satellites;i++){                            // loop over rows (number of satellites in file)
               for (j=0;j<number_of_coords;j++){                            // loop over columns/coordinates in space
                   fscanf(virt_sat_file,"%E",&temp_number);                 // read row element
                   current_sat_coordinates[j] = temp_number;
                   if ((k=fgetc(virt_sat_file))=='\n'){                     // if new line then goto next row (set of spatial coordinates)
                      break;
                   }	
               }
               // vector difference
               current_vector_difference[0] = desired_sat_coordinates[0] - current_sat_coordinates[0];
               current_vector_difference[1] = desired_sat_coordinates[1] - current_sat_coordinates[1];
               // modulus of vector difference
               mod_current_vector_difference = sqrt( current_vector_difference[0] * current_vector_difference[0]
                                                   + current_vector_difference[1] * current_vector_difference[1] );
               mod_closest_vector_difference = sqrt( closest_vector_difference[0] * closest_vector_difference[0]
                                                   + closest_vector_difference[1] * closest_vector_difference[1] );
               // compare the vectors and set the closest satellite vector position
               if (mod_current_vector_difference <= mod_closest_vector_difference){
                  closest_vector_difference[0] = current_vector_difference[0];
                  closest_vector_difference[1] = current_vector_difference[1];
                  closest_sat_coordinates[0] = current_sat_coordinates[0];
                  closest_sat_coordinates[1] = current_sat_coordinates[1];
                  satellite_number_best_match = i+1;                        // +1 because counting starts from zero
                  VirtualSatelliteFilenameFoundString = temp_string;
               }
           }
           fclose(virt_sat_file);                                           // close the current file
        }
    }

    // count the numeber of lines in 'VirtualSatelliteFilenameFoundString' and take only the number of data rows in it
    FILE *file_to_count_lines = fopen(VirtualSatelliteFilenameFoundString.c_str(),"r");
    int line_count=0; while(!fscanf(file_to_count_lines,"%*[^\n]%*c"))line_count++;fseek(file_to_count_lines,0,SEEK_SET);
    fclose(file_to_count_lines);
    line_count = line_count - number_of_satellites;

    // print some messages on screen
    cout << "---> Found satellite No " << satellite_number_best_match << " in file " << "'" << VirtualSatelliteFilenameFoundString << "'" << endl;
    cout << "---> Desired Sat. Coordinates " << desired_sat_coordinates[0] << " " << desired_sat_coordinates[1] << " " << desired_sat_coordinates[2] << endl ;
    cout << "---> Closest Sat. Coordinates " << closest_sat_coordinates[0] << " " << closest_sat_coordinates[1] << " " << closest_sat_coordinates[2] << endl ;
    cout << "---> Closest Vector Difference " << closest_vector_difference[0] << " " << closest_vector_difference[1] << " " << closest_vector_difference[2] << endl ;
    cout << "---> Number of data lines in the file " << line_count << endl;

    // write out best satellite info to be used by read-virt-sat.cpp
    FILE *best_satellite_match_info_found = fopen("BestVirtualSatelliteInfoFoundPoint.txt","w");
    fprintf(best_satellite_match_info_found,"%s %i %i",VirtualSatelliteFilenameFoundString.c_str(),satellite_number_best_match,line_count);
    fclose(best_satellite_match_info_found);

}

////////////////////////////////////////////////////////////////////////
// Function 'read_virtual_sat' uses output from function 'find_virtual_sat' to extract the corresponding satellite data
////////////////////////////////////////////////////////////////////////

void read_virtual_sat(float Bx0,
                      float dt,
                      char output_filename[]){

    const int number_of_coords = 2;            // dimension
    const int number_of_satellites = 16;       // satellites per file
    const int number_of_input_variables = 14;  // these are written in the virt-sat files

    // read the name of the VirtSat filename containing the best satellite match and the respective satellite number
    char VirtualSatelliteFilenameFound[] = "";
    int satellite_number_best_match;
    int line_count;
    char info_fname[] = "BestVirtualSatelliteInfoFoundPoint.txt";
    FILE *satellite_info_file;
    satellite_info_file = fopen(info_fname,"r");
    fscanf(satellite_info_file,"%s %i %i",VirtualSatelliteFilenameFound,&satellite_number_best_match,&line_count);
    fclose(satellite_info_file);

    cout << "---> Extracting data for satellite No " << satellite_number_best_match << " from file " << "'" << VirtualSatelliteFilenameFound << "'" << " containing " << line_count << " lines with data" << endl;

    // the input file form which we will read data for satellite number 'satellite_number_best_match'
    const int satellite_number_in_file = satellite_number_best_match;   // satellite number for which we extract data from input file
    const int number_of_data_rows = line_count;                         // define the number of rows (time records)

    // misc definitions
    const int number_of_elements_per_row = number_of_satellites * number_of_input_variables;
//    float temp_array[number_of_data_rows][number_of_elements_per_row];
//  next line is the 1D equivalent of the above 2D array - memory issues
    float *temp_array = (float*)malloc(number_of_data_rows*number_of_elements_per_row * sizeof(float));
    float temp_array_element;
    int i,j,k;

    // reading the data

    // open the input file
    FILE *virt_sat_file;
    virt_sat_file = fopen(VirtualSatelliteFilenameFound,"r");

    // skip the first 'number_of_satellites' rows containing the (x,y,z) coordinates of the satellites
    for (i=0;i<number_of_satellites;i++){                            // loop row
        for (j=0;j<number_of_coords;j++){                            // loop over columns/coordinates in space
            fscanf(virt_sat_file,"%f",&temp_array_element);          // read row element
            if ((k=fgetc(virt_sat_file))=='\n'){                     // if new line then goto next row
               break;
            }	
        }
    }

    // read the file - read row, then next one, etc...
    for (i=0;i<number_of_data_rows;i++){                             // loop row
        for (j=0;j<number_of_elements_per_row;j++){                  // loop column
            fscanf(virt_sat_file,"%E",&temp_array_element);          // read row element
//            temp_array[i][j]=temp_array_element;                     // store it in temp_array
            temp_array[number_of_elements_per_row * (i - 0) + j]=temp_array_element;
            if ((k=fgetc(virt_sat_file))=='\n'){                     // if new line then goto next row
               break;
            }	
        }
    }

    // close the input file
    fclose(virt_sat_file);

    // define quantities which will be written out
    int ielement;

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

    float ax_tot_e[number_of_data_rows];
    float ay_tot_e[number_of_data_rows];
    float az_tot_e[number_of_data_rows];
    float ax_tot_i[number_of_data_rows];
    float ay_tot_i[number_of_data_rows];
    float az_tot_i[number_of_data_rows];

    float a_tot_e[number_of_data_rows];
    float a_tot_i[number_of_data_rows];

    // definitions for min values
    float min_Bx = 1.0e30;
    float min_By = 1.0e30;
    float min_Bz = 1.0e30;
    float min_Ex = 1.0e30;
    float min_Ey = 1.0e30;
    float min_Ez = 1.0e30;
    float min_Jx_tot_e = 1.0e30;
    float min_Jy_tot_e = 1.0e30;
    float min_Jz_tot_e = 1.0e30;
    float min_Jx_tot_i = 1.0e30;
    float min_Jy_tot_i = 1.0e30;
    float min_Jz_tot_i = 1.0e30;
    float min_rho_tot_e = 1.0e30;
    float min_rho_tot_i = 1.0e30;
    float min_Vx_tot_e = 1.0e30;
    float min_Vy_tot_e = 1.0e30;
    float min_Vz_tot_e = 1.0e30;
    float min_Vx_tot_i = 1.0e30;
    float min_Vy_tot_i = 1.0e30;
    float min_Vz_tot_i = 1.0e30;
    float min_E = 1.0e30;
    float min_B = 1.0e30;
    float min_V_tot_e = 1.0e30;
    float min_V_tot_i = 1.0e30;

    float min_ax_tot_e = 1.0e30;
    float min_ay_tot_e = 1.0e30;
    float min_az_tot_e = 1.0e30;
    float min_ax_tot_i = 1.0e30;
    float min_ay_tot_i = 1.0e30;
    float min_az_tot_i = 1.0e30;
    float min_a_tot_e = 1.0e30;
    float min_a_tot_i = 1.0e30;

    // definitions for max values
    float max_Bx = -1.0e30;
    float max_By = -1.0e30;
    float max_Bz = -1.0e30;
    float max_Ex = -1.0e30;
    float max_Ey = -1.0e30;
    float max_Ez = -1.0e30;
    float max_Jx_tot_e = -1.0e30;
    float max_Jy_tot_e = -1.0e30;
    float max_Jz_tot_e = -1.0e30;
    float max_Jx_tot_i = -1.0e30;
    float max_Jy_tot_i = -1.0e30;
    float max_Jz_tot_i = -1.0e30;
    float max_rho_tot_e = -1.0e30;
    float max_rho_tot_i = -1.0e30;
    float max_Vx_tot_e = -1.0e30;
    float max_Vy_tot_e = -1.0e30;
    float max_Vz_tot_e = -1.0e30;
    float max_Vx_tot_i = -1.0e30;
    float max_Vy_tot_i = -1.0e30;
    float max_Vz_tot_i = -1.0e30;
    float max_E = -1.0e30;
    float max_B = -1.0e30;
    float max_V_tot_e = -1.0e30;
    float max_V_tot_i = -1.0e30;

    float max_ax_tot_e = -1.0e30;
    float max_ay_tot_e = -1.0e30;
    float max_az_tot_e = -1.0e30;
    float max_ax_tot_i = -1.0e30;
    float max_ay_tot_i = -1.0e30;
    float max_az_tot_i = -1.0e30;
    float max_a_tot_e = -1.0e30;
    float max_a_tot_i = -1.0e30;

    int i1=0;

    for (i=0;i<number_of_data_rows;i++){                                // loop time records
        ielement = (satellite_number_in_file - 1) * number_of_input_variables;
        // select from each row the 14 input vars corresponding to the selected satellite
        Bx[i]=temp_array[number_of_elements_per_row * (i - 0) + ielement];
        ielement = ielement + 1;
        By[i]=temp_array[number_of_elements_per_row * (i - 0) + ielement];
        ielement = ielement + 1;
        Bz[i]=temp_array[number_of_elements_per_row * (i - 0) + ielement];
        ielement = ielement + 1;
        Ex[i]=temp_array[number_of_elements_per_row * (i - 0) + ielement];
        ielement = ielement + 1;
        Ey[i]=temp_array[number_of_elements_per_row * (i - 0) + ielement];
        ielement = ielement + 1;
        Ez[i]=temp_array[number_of_elements_per_row * (i - 0) + ielement];
        ielement = ielement + 1;
        Jx_tot_e[i]=temp_array[number_of_elements_per_row * (i - 0) + ielement];
        ielement = ielement + 1;
        Jy_tot_e[i]=temp_array[number_of_elements_per_row * (i - 0) + ielement];
        ielement = ielement + 1;
        Jz_tot_e[i]=temp_array[number_of_elements_per_row * (i - 0) + ielement];
        ielement = ielement + 1;
        Jx_tot_i[i]=temp_array[number_of_elements_per_row * (i - 0) + ielement];
        ielement = ielement + 1;
        Jy_tot_i[i]=temp_array[number_of_elements_per_row * (i - 0) + ielement];
        ielement = ielement + 1;
        Jz_tot_i[i]=temp_array[number_of_elements_per_row * (i - 0) + ielement];
        ielement = ielement + 1;
        rho_tot_e[i]=temp_array[number_of_elements_per_row * (i - 0) + ielement];
        ielement = ielement + 1;
        rho_tot_i[i]=temp_array[number_of_elements_per_row * (i - 0) + ielement];
        // compute additional quantities
        Vx_tot_e[i] = Jx_tot_e[i] / rho_tot_e[i];
        Vy_tot_e[i] = Jy_tot_e[i] / rho_tot_e[i];
        Vz_tot_e[i] = Jz_tot_e[i] / rho_tot_e[i];
        Vx_tot_i[i] = Jx_tot_i[i] / rho_tot_i[i];
        Vy_tot_i[i] = Jy_tot_i[i] / rho_tot_i[i];
        Vz_tot_i[i] = Jz_tot_i[i] / rho_tot_i[i];
        B[i] = sqrt( Bx[i]*Bx[i] + By[i]*By[i] + Bz[i]*Bz[i] );
        E[i] = sqrt( Ex[i]*Ex[i] + Ey[i]*Ey[i] + Ez[i]*Ez[i] );
        V_tot_e[i]=sqrt(Vx_tot_e[i]*Vx_tot_e[i]+Vy_tot_e[i]*Vy_tot_e[i]+Vz_tot_e[i]*Vz_tot_e[i]);
        V_tot_i[i]=sqrt(Vx_tot_i[i]*Vx_tot_i[i]+Vy_tot_i[i]*Vy_tot_i[i]+Vz_tot_i[i]*Vz_tot_i[i]);
        // accelerations
        if (i>0) i1=i-1;
        ax_tot_e[i] = (Vx_tot_e[i]-Vx_tot_e[i1]) / dt;
        ay_tot_e[i] = (Vy_tot_e[i]-Vy_tot_e[i1]) / dt;
        az_tot_e[i] = (Vz_tot_e[i]-Vz_tot_e[i1]) / dt;
        ax_tot_i[i] = (Vx_tot_i[i]-Vx_tot_i[i1]) / dt;
        ay_tot_i[i] = (Vy_tot_i[i]-Vy_tot_i[i1]) / dt;
        az_tot_i[i] = (Vz_tot_i[i]-Vz_tot_i[i1]) / dt;
        a_tot_e[i]=sqrt(ax_tot_e[i]*ax_tot_e[i]+ay_tot_e[i]*ay_tot_e[i]+az_tot_e[i]*az_tot_e[i]);
        a_tot_i[i]=sqrt(ax_tot_i[i]*ax_tot_i[i]+ay_tot_i[i]*ay_tot_i[i]+az_tot_i[i]*az_tot_i[i]);

        // check for min and max for each quantity
        if (Bx[i]<=min_Bx) {
           min_Bx=Bx[i];
        }
        if (Bx[i]>=max_Bx) {
           max_Bx=Bx[i];
        }
        if (By[i]<=min_By) {
           min_By=By[i];
        }
        if (By[i]>=max_By) {
           max_By=By[i];
        }
        if (Bz[i]<=min_Bz) {
           min_Bz=Bz[i];
        }
        if (Bz[i]>=max_Bz) {
           max_Bz=Bz[i];
        }
        if (Ex[i]<=min_Ex) {
           min_Ex=Ex[i];
        }
        if (Ex[i]>=max_Ex) {
           max_Ex=Ex[i];
        }
        if (Ey[i]<=min_Ey) {
           min_Ey=Ey[i];
        }
        if (Ey[i]>=max_Ey) {
           max_Ey=Ey[i];
        }
        if (Ez[i]<=min_Ez) {
           min_Ez=Ez[i];
        }
        if (Ez[i]>=max_Ez) {
           max_Ez=Ez[i];
        }
        if (Vx_tot_e[i]<=min_Vx_tot_e) {
           min_Vx_tot_e=Vx_tot_e[i];
        }
        if (Vx_tot_e[i]>=max_Vx_tot_e) {
           max_Vx_tot_e=Vx_tot_e[i];
        }
        if (Vy_tot_e[i]<=min_Vy_tot_e) {
           min_Vy_tot_e=Vy_tot_e[i];
        }
        if (Vy_tot_e[i]>=max_Vy_tot_e) {
           max_Vy_tot_e=Vy_tot_e[i];
        }
        if (Vz_tot_e[i]<=min_Vz_tot_e) {
           min_Vz_tot_e=Vz_tot_e[i];
        }
        if (Vz_tot_e[i]>=max_Vz_tot_e) {
           max_Vz_tot_e=Vz_tot_e[i];
        }
        if (Vx_tot_i[i]<=min_Vx_tot_i) {
           min_Vx_tot_i=Vx_tot_i[i];
        }
        if (Vx_tot_i[i]>=max_Vx_tot_i) {
           max_Vx_tot_i=Vx_tot_i[i];
        }
        if (Vy_tot_i[i]<=min_Vy_tot_i) {
           min_Vy_tot_i=Vy_tot_i[i];
        }
        if (Vy_tot_i[i]>=max_Vy_tot_i) {
           max_Vy_tot_i=Vy_tot_i[i];
        }
        if (Vz_tot_i[i]<=min_Vz_tot_i) {
           min_Vz_tot_i=Vz_tot_i[i];
        }
        if (Vz_tot_i[i]>=max_Vz_tot_i) {
           max_Vz_tot_i=Vz_tot_i[i];
        }
        if (Jx_tot_e[i]<=min_Jx_tot_e) {
           min_Jx_tot_e=Jx_tot_e[i];
        }
        if (Jx_tot_e[i]>=max_Jx_tot_e) {
           max_Jx_tot_e=Jx_tot_e[i];
        }
        if (Jy_tot_e[i]<=min_Jy_tot_e) {
           min_Jy_tot_e=Jy_tot_e[i];
        }
        if (Jy_tot_e[i]>=max_Jy_tot_e) {
           max_Jy_tot_e=Jy_tot_e[i];
        }
        if (Jz_tot_e[i]<=min_Jz_tot_e) {
           min_Jz_tot_e=Jz_tot_e[i];
        }
        if (Jz_tot_e[i]>=max_Jz_tot_e) {
           max_Jz_tot_e=Jz_tot_e[i];
        }
        if (Jx_tot_i[i]<=min_Jx_tot_i) {
           min_Jx_tot_i=Jx_tot_i[i];
        }
        if (Jx_tot_i[i]>=max_Jx_tot_i) {
           max_Jx_tot_i=Jx_tot_i[i];
        }
        if (Jy_tot_i[i]<=min_Jy_tot_i) {
           min_Jy_tot_i=Jy_tot_i[i];
        }
        if (Jy_tot_i[i]>=max_Jy_tot_i) {
           max_Jy_tot_i=Jy_tot_i[i];
        }
        if (Jz_tot_i[i]<=min_Jz_tot_i) {
           min_Jz_tot_i=Jz_tot_i[i];
        }
        if (Jz_tot_i[i]>=max_Jz_tot_i) {
           max_Jz_tot_i=Jz_tot_i[i];
        }
        if (rho_tot_e[i]<=min_rho_tot_e) {
           min_rho_tot_e=rho_tot_e[i];
        }
        if (rho_tot_e[i]>=max_rho_tot_e) {
           max_rho_tot_e=rho_tot_e[i];
        }
        if (rho_tot_i[i]<=min_rho_tot_i) {
           min_rho_tot_i=rho_tot_i[i];
        }
        if (rho_tot_i[i]>=max_rho_tot_i) {
           max_rho_tot_i=rho_tot_i[i];
        }
        if (V_tot_e[i]<=min_V_tot_e) {
           min_V_tot_e=V_tot_e[i];
        }
        if (V_tot_e[i]>=max_V_tot_e) {
           max_V_tot_e=V_tot_e[i];
        }
        if (V_tot_i[i]<=min_V_tot_i) {
           min_V_tot_i=V_tot_i[i];
        }
        if (V_tot_i[i]>=max_V_tot_i) {
           max_V_tot_i=V_tot_i[i];
        }
        if (B[i]<=min_B) {
           min_B=B[i];
        }
        if (B[i]>=max_B) {
           max_B=B[i];
        }
        if (E[i]<=min_E) {
           min_E=E[i];
        }
        if (E[i]>=max_E) {
           max_E=E[i];
        }
        // accelerations
        if (ax_tot_e[i]<=min_ax_tot_e) {
           min_ax_tot_e=ax_tot_e[i];
        }
        if (ax_tot_e[i]>=max_ax_tot_e) {
           max_ax_tot_e=ax_tot_e[i];
        }
        if (ay_tot_e[i]<=min_ay_tot_e) {
           min_ay_tot_e=ay_tot_e[i];
        }
        if (ay_tot_e[i]>=max_ay_tot_e) {
           max_ay_tot_e=ay_tot_e[i];
        }
        if (az_tot_e[i]<=min_az_tot_e) {
           min_az_tot_e=az_tot_e[i];
        }
        if (az_tot_e[i]>=max_az_tot_e) {
           max_az_tot_e=az_tot_e[i];
        }
        if (ax_tot_i[i]<=min_ax_tot_i) {
           min_ax_tot_i=ax_tot_i[i];
        }
        if (ax_tot_i[i]>=max_ax_tot_i) {
           max_ax_tot_i=ax_tot_i[i];
        }
        if (ay_tot_i[i]<=min_ay_tot_i) {
           min_ay_tot_i=ay_tot_i[i];
        }
        if (ay_tot_i[i]>=max_ay_tot_i) {
           max_ay_tot_i=ay_tot_i[i];
        }
        if (az_tot_i[i]<=min_az_tot_i) {
           min_az_tot_i=az_tot_i[i];
        }
        if (az_tot_i[i]>=max_az_tot_i) {
           max_az_tot_i=az_tot_i[i];
        }
        if (a_tot_e[i]<=min_a_tot_e) {
           min_a_tot_e=a_tot_e[i];
        }
        if (a_tot_e[i]>=max_a_tot_e) {
           max_a_tot_e=a_tot_e[i];
        }
        if (a_tot_i[i]<=min_a_tot_i) {
           min_a_tot_i=a_tot_i[i];
        }
        if (a_tot_i[i]>=max_a_tot_i) {
           max_a_tot_i=a_tot_i[i];
        }

    }

    delete[] temp_array;

    // writing the data in file 'output_filename'

    // open the output file
    FILE *output_file;
    output_file = fopen(output_filename,"w");

    float omega_ci_x_t;                                                 // this will be the X-axis in the output file

    // write out notation at the beginning
    fprintf(output_file,"%s","omega_ci_x_t ");
    fprintf(output_file,"%s","Bx            By            Bz            Ex            Ey            Ez");
    fprintf(output_file,"%s","     Jx_tot_e      Jy_tot_e      Jz_tot_e      Jx_tot_i");
    fprintf(output_file,"%s","      Jy_tot_i      Jz_tot_i    rho_tot_e    rho_tot_i      ");
    fprintf(output_file,"%s","Vx_tot_e     Vy_tot_e     Vz_tot_e     Vx_tot_i     Vy_tot_i     Vz_tot_i");
    fprintf(output_file,"%s","          B            E            ");
    fprintf(output_file,"%s","   V_tot_e     V_tot_i");
    fprintf(output_file,"%s","ax_tot_e     ay_tot_e     az_tot_e     ax_tot_i     ay_tot_i     az_tot_i");
    fprintf(output_file,"%s","   a_tot_e     a_tot_i");
    fprintf(output_file,"\n","");

    // write out the data
    for (i=0;i<number_of_data_rows;i++){                                // loop time records
        omega_ci_x_t = Bx0*i*dt;
        fprintf(output_file,"%f ",omega_ci_x_t);
        fprintf(output_file,"%E %E %E %E %E %E %E %E %E %E %E %E %E %E "
                           ,Bx[i],By[i],Bz[i]
                           ,Ex[i],Ey[i],Ez[i]
                           ,Jx_tot_e[i],Jy_tot_e[i],Jz_tot_e[i]
                           ,Jx_tot_i[i],Jy_tot_i[i],Jz_tot_i[i]
                           ,rho_tot_e[i],rho_tot_i[i]);
        fprintf(output_file,"%E %E %E %E %E %E "
                           ,Vx_tot_e[i],Vy_tot_e[i],Vz_tot_e[i]
                           ,Vx_tot_i[i],Vy_tot_i[i],Vz_tot_i[i]);
        fprintf(output_file,"%E %E "
                           ,B[i],E[i]);
        fprintf(output_file,"%E %E "
                           ,V_tot_e[i],V_tot_i[i]);
        fprintf(output_file,"%E %E %E %E %E %E "
                           ,ax_tot_e[i],ay_tot_e[i],az_tot_e[i]
                           ,ax_tot_i[i],ay_tot_i[i],az_tot_i[i]);
        fprintf(output_file,"%E %E "
                           ,a_tot_e[i],a_tot_i[i]);
        fprintf(output_file,"\n","");
    }

    // close the output file
    fclose(output_file);

    // create minmax files, write out min and max, close the files
    // create
    FILE *output_file_Bx_minmax = fopen("output_filename_Bx_minmax.txt","w");
    FILE *output_file_By_minmax = fopen("output_filename_By_minmax.txt","w");
    FILE *output_file_Bz_minmax = fopen("output_filename_Bz_minmax.txt","w");
    FILE *output_file_Ex_minmax = fopen("output_filename_Ex_minmax.txt","w");
    FILE *output_file_Ey_minmax = fopen("output_filename_Ey_minmax.txt","w");
    FILE *output_file_Ez_minmax = fopen("output_filename_Ez_minmax.txt","w");
    FILE *output_file_Jx_tot_e_minmax = fopen("output_filename_Jx_tot_e_minmax.txt","w");
    FILE *output_file_Jy_tot_e_minmax = fopen("output_filename_Jy_tot_e_minmax.txt","w");
    FILE *output_file_Jz_tot_e_minmax = fopen("output_filename_Jz_tot_e_minmax.txt","w");
    FILE *output_file_Jx_tot_i_minmax = fopen("output_filename_Jx_tot_i_minmax.txt","w");
    FILE *output_file_Jy_tot_i_minmax = fopen("output_filename_Jy_tot_i_minmax.txt","w");
    FILE *output_file_Jz_tot_i_minmax = fopen("output_filename_Jz_tot_i_minmax.txt","w");
    FILE *output_file_rho_tot_e_minmax = fopen("output_filename_rho_tot_e_minmax.txt","w");
    FILE *output_file_rho_tot_i_minmax = fopen("output_filename_rho_tot_i_minmax.txt","w");
    FILE *output_file_Vx_tot_e_minmax = fopen("output_filename_Vx_tot_e_minmax.txt","w");
    FILE *output_file_Vy_tot_e_minmax = fopen("output_filename_Vy_tot_e_minmax.txt","w");
    FILE *output_file_Vz_tot_e_minmax = fopen("output_filename_Vz_tot_e_minmax.txt","w");
    FILE *output_file_Vx_tot_i_minmax = fopen("output_filename_Vx_tot_i_minmax.txt","w");
    FILE *output_file_Vy_tot_i_minmax = fopen("output_filename_Vy_tot_i_minmax.txt","w");
    FILE *output_file_Vz_tot_i_minmax = fopen("output_filename_Vz_tot_i_minmax.txt","w");
    FILE *output_file_B_minmax = fopen("output_filename_B_minmax.txt","w");
    FILE *output_file_E_minmax = fopen("output_filename_E_minmax.txt","w");
    FILE *output_file_V_tot_e_minmax = fopen("output_filename_V_tot_e_minmax.txt","w");
    FILE *output_file_V_tot_i_minmax = fopen("output_filename_V_tot_i_minmax.txt","w");
    FILE *output_file_ax_tot_e_minmax = fopen("output_filename_ax_tot_e_minmax.txt","w");
    FILE *output_file_ay_tot_e_minmax = fopen("output_filename_ay_tot_e_minmax.txt","w");
    FILE *output_file_az_tot_e_minmax = fopen("output_filename_az_tot_e_minmax.txt","w");
    FILE *output_file_ax_tot_i_minmax = fopen("output_filename_ax_tot_i_minmax.txt","w");
    FILE *output_file_ay_tot_i_minmax = fopen("output_filename_ay_tot_i_minmax.txt","w");
    FILE *output_file_az_tot_i_minmax = fopen("output_filename_az_tot_i_minmax.txt","w");
    FILE *output_file_a_tot_e_minmax = fopen("output_filename_a_tot_e_minmax.txt","w");
    FILE *output_file_a_tot_i_minmax = fopen("output_filename_a_tot_i_minmax.txt","w");
    // write
    fprintf(output_file_Bx_minmax,"%E %E ",min_Bx,max_Bx);
    fprintf(output_file_By_minmax,"%E %E ",min_By,max_By);
    fprintf(output_file_Bz_minmax,"%E %E ",min_Bz,max_Bz);
    fprintf(output_file_Ex_minmax,"%E %E ",min_Ex,max_Ex);
    fprintf(output_file_Ey_minmax,"%E %E ",min_Ey,max_Ey);
    fprintf(output_file_Ez_minmax,"%E %E ",min_Ez,max_Ez);
    fprintf(output_file_Jx_tot_e_minmax,"%E %E ",min_Jx_tot_e,max_Jx_tot_e);
    fprintf(output_file_Jy_tot_e_minmax,"%E %E ",min_Jy_tot_e,max_Jy_tot_e);
    fprintf(output_file_Jz_tot_e_minmax,"%E %E ",min_Jz_tot_e,max_Jz_tot_e);
    fprintf(output_file_Jx_tot_i_minmax,"%E %E ",min_Jx_tot_i,max_Jx_tot_i);
    fprintf(output_file_Jy_tot_i_minmax,"%E %E ",min_Jy_tot_i,max_Jy_tot_i);
    fprintf(output_file_Jz_tot_i_minmax,"%E %E ",min_Jz_tot_i,max_Jz_tot_i);
    fprintf(output_file_rho_tot_e_minmax,"%E %E ",min_rho_tot_e,max_rho_tot_e);
    fprintf(output_file_rho_tot_i_minmax,"%E %E ",min_rho_tot_i,max_rho_tot_i);
    fprintf(output_file_Vx_tot_e_minmax,"%E %E ",min_Vx_tot_e,max_Vx_tot_e);
    fprintf(output_file_Vy_tot_e_minmax,"%E %E ",min_Vy_tot_e,max_Vy_tot_e);
    fprintf(output_file_Vz_tot_e_minmax,"%E %E ",min_Vz_tot_e,max_Vz_tot_e);
    fprintf(output_file_Vx_tot_i_minmax,"%E %E ",min_Vx_tot_i,max_Vx_tot_i);
    fprintf(output_file_Vy_tot_i_minmax,"%E %E ",min_Vy_tot_i,max_Vy_tot_i);
    fprintf(output_file_Vz_tot_i_minmax,"%E %E ",min_Vz_tot_i,max_Vz_tot_i);
    fprintf(output_file_B_minmax,"%E %E ",min_B,max_B);
    fprintf(output_file_E_minmax,"%E %E ",min_E,max_E);
    fprintf(output_file_V_tot_e_minmax,"%E %E ",min_V_tot_e,max_V_tot_e);
    fprintf(output_file_V_tot_i_minmax,"%E %E ",min_V_tot_i,max_V_tot_i);
    fprintf(output_file_ax_tot_e_minmax,"%E %E ",min_ax_tot_e,max_ax_tot_e);
    fprintf(output_file_ay_tot_e_minmax,"%E %E ",min_ay_tot_e,max_ay_tot_e);
    fprintf(output_file_az_tot_e_minmax,"%E %E ",min_az_tot_e,max_az_tot_e);
    fprintf(output_file_ax_tot_i_minmax,"%E %E ",min_ax_tot_i,max_ax_tot_i);
    fprintf(output_file_ay_tot_i_minmax,"%E %E ",min_ay_tot_i,max_ay_tot_i);
    fprintf(output_file_az_tot_i_minmax,"%E %E ",min_az_tot_i,max_az_tot_i);
    fprintf(output_file_a_tot_e_minmax,"%E %E ",min_a_tot_e,max_a_tot_e);
    fprintf(output_file_a_tot_i_minmax,"%E %E ",min_a_tot_i,max_a_tot_i);
    // close
    fclose(output_file_Bx_minmax);
    fclose(output_file_By_minmax);
    fclose(output_file_Bz_minmax);
    fclose(output_file_Ex_minmax);
    fclose(output_file_Ey_minmax);
    fclose(output_file_Ez_minmax);
    fclose(output_file_Jx_tot_e_minmax);
    fclose(output_file_Jy_tot_e_minmax);
    fclose(output_file_Jz_tot_e_minmax);
    fclose(output_file_Jx_tot_i_minmax);
    fclose(output_file_Jy_tot_i_minmax);
    fclose(output_file_Jz_tot_i_minmax);
    fclose(output_file_rho_tot_e_minmax);
    fclose(output_file_rho_tot_i_minmax);
    fclose(output_file_Vx_tot_e_minmax);
    fclose(output_file_Vy_tot_e_minmax);
    fclose(output_file_Vz_tot_e_minmax);
    fclose(output_file_Vx_tot_i_minmax);
    fclose(output_file_Vy_tot_i_minmax);
    fclose(output_file_Vz_tot_i_minmax);
    fclose(output_file_B_minmax);
    fclose(output_file_E_minmax);
    fclose(output_file_V_tot_e_minmax);
    fclose(output_file_V_tot_i_minmax);
    fclose(output_file_ax_tot_e_minmax);
    fclose(output_file_ay_tot_e_minmax);
    fclose(output_file_az_tot_e_minmax);
    fclose(output_file_ax_tot_i_minmax);
    fclose(output_file_ay_tot_i_minmax);
    fclose(output_file_az_tot_i_minmax);
    fclose(output_file_a_tot_e_minmax);
    fclose(output_file_a_tot_i_minmax);

    cout << "---> Finished data extraction." << endl;

}
// END FUNCTION DEFINITION

// BEGIN MAIN
int main()
{

    // user inputs
    #include "user_inputs_point_2D.h"
    // find the virtual satellite which is closest to our specified (x,y,z) point
    find_virtual_sat (x_satellite_user,
                      y_satellite_user,
                      nproc);
    // extract the data for that virtual satellite
    read_virtual_sat (Bx0,dt,output_filename);

    return 0;

}
// END MAIN


