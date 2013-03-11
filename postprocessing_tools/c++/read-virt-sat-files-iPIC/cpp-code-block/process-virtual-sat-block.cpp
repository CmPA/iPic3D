/***************************************************************************

    process-virtual-sat.cpp  -  postprocessing program to extract data from VirtualSatellites* files from iPIC3D

    Author:  A.E.Vapirev, KU Leuven, Afdeling Plasma-astrofysica
             2011 Sep 3

    Change log           :

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
// Function 'find_virtual_sat_block' finds the best satellite match with closest position to user given (x,y,z) spatial coordinates.
// The result is used by funciton 'read_virtual_sat_block'.
////////////////////////////////////////////////////////////////////////

void find_virtual_sat_block(float x_satellite_user,
                            float y_satellite_user,
                            float z_satellite_user,
                            int const nproc,
                            string XYZ_block){

    FILE *virt_sat_file;                       // input file

    const int number_of_coords = 3;            // dimension
    const int number_of_satellites = 27;       // satellites per file

    int i,j,k,iproc;
    float temp_number[number_of_coords];

    int satellite_number_best_match[nproc];                                   // later used in read-virtual-sat.cpp to extract the corresponding data
    string VirtualSatelliteFilenameFoundString[nproc];                        // this file contains the best satellite coordinates match
    int iSatCount = 0;

    vector<float> desired_sat_coordinates (number_of_coords-1,0);             // define 3D vector with zeros - the user desired satellite position
    vector<float> current_sat_coordinates (number_of_coords-1,0);             // define 3D vector with zeros - reading current coords. from the file
    vector<float> closest_sat_coordinates (number_of_coords-1,1.0E10);        // define 3D vector - closest satellite position to what we want
    vector<float> current_vector_difference (number_of_coords-1,0); 
    vector<float> closest_vector_difference (number_of_coords-1,1.0E10);      // this will end up being closest to the desired satellite position

    float mod_current_vector_difference;
    float mod_closest_vector_difference;

    if ( XYZ_block == "XY" ){
       // set up vector with the user defined coordinates
       desired_sat_coordinates[0] = x_satellite_user;
       desired_sat_coordinates[1] = y_satellite_user;
    }
    if ( XYZ_block == "XZ" ){
       // set up vector with the user defined coordinates
       desired_sat_coordinates[0] = x_satellite_user;
       desired_sat_coordinates[1] = z_satellite_user;
    }
    if ( XYZ_block == "YZ" ){
       // set up vector with the user defined coordinates
       desired_sat_coordinates[0] = y_satellite_user;
       desired_sat_coordinates[1] = z_satellite_user;
    }

    float original_desired_sat_coordinates0=desired_sat_coordinates[0];
    float original_desired_sat_coordinates1=desired_sat_coordinates[1];

    string temp_string;

    cout << "---> Searching for best satellite match..." << endl;

    // start best satellite search over all files
    for (int iproc=0; iproc < nproc; iproc++){                              // loop over all files
        stringstream ss;
        ss << iproc;
        temp_string = "VirtualSatelliteTraces" + ss.str() + ".txt";
        virt_sat_file = fopen(temp_string.c_str(),"r");

        // first scan for the closest coordinates to the desired ones and store them
        if (virt_sat_file){                                                 // if the file exists
           for (i=0;i<number_of_satellites;i++){                            // loop over rows (number of satellites in file)
               for (j=0;j<number_of_coords;j++){                            // loop over columns/coordinates in space
                   fscanf(virt_sat_file,"%E",&temp_number[j]);              // read row element
                   if ((k=fgetc(virt_sat_file))=='\n'){                     // if new line then goto next row (set of spatial coordinates)
                      break;
                   }
               }
               // set the respective current_sat_coordinates
               if ( XYZ_block == "XY" ){
                  current_sat_coordinates[0] = temp_number[0];
                  current_sat_coordinates[1] = temp_number[1];
               }
               if ( XYZ_block == "XZ" ){
                  current_sat_coordinates[0] = temp_number[0];
                  current_sat_coordinates[1] = temp_number[2];
               }
               if ( XYZ_block == "YZ" ){
                  current_sat_coordinates[0] = temp_number[1];
                  current_sat_coordinates[1] = temp_number[2];
               }
               // vector difference
               current_vector_difference[0] = desired_sat_coordinates[0] - current_sat_coordinates[0];
               current_vector_difference[1] = desired_sat_coordinates[1] - current_sat_coordinates[1];
               // modulus of vector difference
               mod_current_vector_difference = sqrt( current_vector_difference[0] * current_vector_difference[0]
                                                   + current_vector_difference[1] * current_vector_difference[1]);
               mod_closest_vector_difference = sqrt( closest_vector_difference[0] * closest_vector_difference[0]
                                                   + closest_vector_difference[1] * closest_vector_difference[1]);
               // check and set up the closest coordinates
               if (mod_current_vector_difference <= mod_closest_vector_difference){
                      closest_vector_difference[0] = current_vector_difference[0];
                      closest_vector_difference[1] = current_vector_difference[1];
                      closest_sat_coordinates[0] = current_sat_coordinates[0];
                      closest_sat_coordinates[1] = current_sat_coordinates[1];
               }
           }
           fclose(virt_sat_file);                                           // close the current file
        }

    }

   desired_sat_coordinates[0]=closest_sat_coordinates[0];
   desired_sat_coordinates[1]=closest_sat_coordinates[1];

    // now we have found the closest coordinates to the desired ones and search through all satellites for that coordinate trace in space
    for (int iproc=0; iproc < nproc; iproc++){                              // loop over all files
        stringstream ss;
        ss << iproc;
        temp_string = "VirtualSatelliteTraces" + ss.str() + ".txt";
        virt_sat_file = fopen(temp_string.c_str(),"r");

        // scan the file again and this time when coordinate match is found, record the satellite number for later data extraction
        if (virt_sat_file){                                                 // if the file exists
           for (i=0;i<number_of_satellites;i++){                            // loop over rows (number of satellites in file)
               for (j=0;j<number_of_coords;j++){                            // loop over columns/coordinates in space
                   fscanf(virt_sat_file,"%E",&temp_number[j]);              // read row element
                   if ((k=fgetc(virt_sat_file))=='\n'){                     // if new line then goto next row (set of spatial coordinates)
                      break;
                   }
               }
               // set the respective current_sat_coordinates
               if ( XYZ_block == "XY" ){
                  current_sat_coordinates[0] = temp_number[0];
                  current_sat_coordinates[1] = temp_number[1];
               }
               if ( XYZ_block == "XZ" ){
                  current_sat_coordinates[0] = temp_number[0];
                  current_sat_coordinates[1] = temp_number[2];
               }
               if ( XYZ_block == "YZ" ){
                  current_sat_coordinates[0] = temp_number[1];
                  current_sat_coordinates[1] = temp_number[2];
               }
               // vector difference
               current_vector_difference[0] = desired_sat_coordinates[0] - current_sat_coordinates[0];
               current_vector_difference[1] = desired_sat_coordinates[1] - current_sat_coordinates[1];
               // modulus of vector difference
               mod_current_vector_difference = sqrt( current_vector_difference[0] * current_vector_difference[0]
                                                   + current_vector_difference[1] * current_vector_difference[1]);
               mod_closest_vector_difference = sqrt( closest_vector_difference[0] * closest_vector_difference[0]
                                                   + closest_vector_difference[1] * closest_vector_difference[1]);
               // compare the vectors and set the satellite number in the file
               if (mod_current_vector_difference <= mod_closest_vector_difference){
                  satellite_number_best_match[iSatCount] = i+1;             // +1 because counting starts from zero
                  VirtualSatelliteFilenameFoundString[iSatCount] = temp_string;
                  iSatCount = iSatCount + 1;
               }
           }
           fclose(virt_sat_file);                                           // close the current file
        }
    }

    // write out best satellite info in a file 'BestVirtualSatelliteInfoFoundBlock.txt' to be used by read-virt-sat.cpp
    for (i=0;i<iSatCount;i++){
        // count the numeber of lines in 'VirtualSatelliteFilenameFoundString[i]' and take only the number of data rows in it
        FILE *file_to_count_lines = fopen(VirtualSatelliteFilenameFoundString[i].c_str(),"r");
        int line_count=0; while(!fscanf(file_to_count_lines,"%*[^\n]%*c"))line_count++;fseek(file_to_count_lines,0,SEEK_SET);
        fclose(file_to_count_lines);
        line_count = line_count - number_of_satellites;

        // print some messages on screen
        cout << endl;
        cout << "---> Found satellite No " << satellite_number_best_match[i] << " in file " << "'" << VirtualSatelliteFilenameFoundString[i] << "'" << endl;
        cout << "---> Desired " << XYZ_block << " Sat. Coordinates " << original_desired_sat_coordinates0 << " " << original_desired_sat_coordinates1 << endl ;
        cout << "---> Closest " << XYZ_block << " Sat. Coordinates " << closest_sat_coordinates[0] << " " << closest_sat_coordinates[1] << endl ;
        cout << "---> Closest " << XYZ_block << " Vector Difference " << closest_vector_difference[0] << " " << closest_vector_difference[1] << endl ;
        cout << "---> Number of data lines in the file " << line_count << endl;

        // write out best satellite info to be used by read-virt-sat.cpp
        FILE *best_satellite_match_info_found = fopen("BestVirtualSatelliteInfoFoundBlock.txt","a");
        fprintf(best_satellite_match_info_found,"%s %i %i",VirtualSatelliteFilenameFoundString[i].c_str(),satellite_number_best_match[i],line_count);
        fprintf(best_satellite_match_info_found,"\n","");
        fclose(best_satellite_match_info_found);
    }

}

////////////////////////////////////////////////////////////////////////
// Function'read_virtual_sat_block' uses output from function 'find_virtual_sat_block' to extract the corresponding satellite data
////////////////////////////////////////////////////////////////////////

void read_virtual_sat_block(float Bx0,
                            float dt,
                            char output_filename[]){

    const int number_of_coords = 3;            // dimension
    const int number_of_satellites = 27;       // satellites per file
    const int number_of_input_variables = 14;  // these are written in the virt-sat files

    // count the numeber of lines in file 'BestVirtualSatelliteInfoFoundBlock.txt'
    FILE *file_to_count_lines = fopen("BestVirtualSatelliteInfoFoundBlock.txt","r");
    int number_of_satellites_to_extract=0; while(!fscanf(file_to_count_lines,"%*[^\n]%*c"))number_of_satellites_to_extract++;fseek(file_to_count_lines,0,SEEK_SET);
    fclose(file_to_count_lines);

    // open the satellite info file
    char VirtualSatelliteFilenameFound[] = "";
    int satellite_number_best_match;
    int line_count;
    char info_fname[] = "BestVirtualSatelliteInfoFoundBlock.txt";
    FILE *satellite_info_file;
    satellite_info_file = fopen(info_fname,"r");

    // open the output file
    FILE *output_file;
    output_file = fopen(output_filename,"w");

    // loop over all satellites from 'BestVirtualSatelliteInfoFoundBlock.txt' and extract the data for each of them
    for (int isatloop=0;isatloop<number_of_satellites_to_extract;isatloop++){

    // read the name of the VirtSat filename containing the best satellite match and the respective satellite number
    fscanf(satellite_info_file,"%s %i %i",VirtualSatelliteFilenameFound,&satellite_number_best_match,&line_count);

    cout << endl;
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

    float dBdt[number_of_data_rows];
    float dEdt[number_of_data_rows];

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
    }

    delete[] temp_array;

    // write out - each block of data corresponds to a satellite 'VirtualSatelliteFilenameFound'

    float omega_ci_x_t;                                                 // this will be the X-axis in the output file

    // write out notation at the beginning of each block
    fprintf(output_file,"%s %i %s %s","Data Block: satellite No",satellite_number_best_match,"from file",VirtualSatelliteFilenameFound);
    fprintf(output_file,"\n","");
    fprintf(output_file,"%i %s",number_of_data_rows,"data rows in this block");
    fprintf(output_file,"\n","");
    fprintf(output_file,"%s","omega_ci_x_t ");
    fprintf(output_file,"%s","Bx            By            Bz            Ex            Ey            Ez");
    fprintf(output_file,"%s","     Jx_tot_e      Jy_tot_e      Jz_tot_e      Jx_tot_i");
    fprintf(output_file,"%s","      Jy_tot_i      Jz_tot_i    rho_tot_e    rho_tot_i      ");
    fprintf(output_file,"%s","Vx_tot_e     Vy_tot_e     Vz_tot_e     Vx_tot_i     Vy_tot_i     Vz_tot_i");
    fprintf(output_file,"%s","          B            E            ");
    fprintf(output_file,"%s","   V_tot_e     V_tot_i");
    fprintf(output_file,"%s","dBdt         dEdt          dBxdt        dBydt        dBzdt        dExdt        dEydt        dEzdt      ");
    fprintf(output_file,"%s","FFTBx         FFTBy        FFTBz         FFTEx         FFTEy         FFTEz");
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
        fprintf(output_file,"\n","");
    }

    cout << "---> Finished data extraction." << endl;

    } // end isatloop over the selected satellites

    fclose(output_file);                                                // close the output file
    fclose(satellite_info_file);                                        // close file 'BestVirtualSatelliteInfoFoundBlock.txt'

}
// END FUNCTION DEFINITION

// BEGIN MAIN
int main()
{

    // user inputs
    #include "user_inputs_block.h"
    // find the trace of virtual satellites which is closest to our specified (XYZ_block) line
    find_virtual_sat_block (x_satellite_user,
                            y_satellite_user,
                            z_satellite_user,
                            nproc,
                            XYZ_block);
    // extract the data for that set of virtual satellites
    read_virtual_sat_block (Bx0,dt,output_filename);

    return 0;

}
// END MAIN


