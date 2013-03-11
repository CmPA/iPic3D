#!/bin/bash
########################################################################
#
#    Author:  A.E.Vapirev, KU Leuven, Afdeling Plasma-astrofysica
#             2011 Sep 3
#
#    Help:
#
#    Change log           :
#
########################################################################

########################## define user inputs ##########################

# inputs for finding the best satellite
x_satellite_user=20.878                            # x-coordinate of virtual satellite
y_satellite_user=4.94141                           # y-coordinate of virtual satellite
nproc=2048                                         # number of files to search into

# inputs to write out the data for the satellite
Bx0=0.003638                                       # base field x
dt=0.060                                           # timestep

# input and output folder
Cpp_code_Folder="/home/alexander/Programs/c/read-virt-sat-files-iPIC/cpp-code-point-2D"
Input_Data_Folder="/home/alexander/data/test/virt-sat-data-folder/datarun157.2"
Output_Data_Folder="/home/alexander/data/test/temp-point/output-data-point-x20-run157.2"

########################### end user inputs ############################

# output file - better to not change it because the filenames in this c++/gnuplot suite are inter-connected
output_filename="sat_output_point.txt"

# check if input folder and output folder are the same. If so, force exit - never a good idea to postprocess in the input directory.
if [ ${Input_Data_Folder} = ${Output_Data_Folder} ]; then
   echo
   echo "---> 'Input_Data_Folder' is the same as 'Output_Data_Folder' in 'user inputs'. It is never a good idea to postprocess in the input directory. Fix the path names and rerun the script. Forcing exit..."
   echo
   exit
fi

# create output folder if it does not exist
if [ -e ${Output_Data_Folder} ]; then
   echo
   echo "---> The 'Output_Data_Folder' exists. Do you want to remove it and create an empty one? - then press 'y' and 'enter'. If you dont want to delete it, press 'n' and 'enter' and then change it with a new one in 'user inputs'."
   echo
   echo "y/n? and then 'enter'"
   echo
   read yes_no
   if [ "${yes_no}" = "y" ]; then
      rm -r ${Output_Data_Folder}
   else
      echo
      echo "---> exiting..."
      exit
   fi
   echo
fi
mkdir ${Output_Data_Folder}

# remove previous user_inputs.h if existing
if [ -e ${Cpp_code_Folder}/user_inputs_point_2D.h ];then
   rm ${Cpp_code_Folder}/user_inputs_point_2D.h
fi

echo

# create user_inputs_2D.h file
echo "    float x_satellite_user = $x_satellite_user;         // x-coordinate"                     >> ${Cpp_code_Folder}/user_inputs_point_2D.h
echo "    float y_satellite_user = $y_satellite_user;         // y-coordinate"                     >> ${Cpp_code_Folder}/user_inputs_point_2D.h
echo "    int const nproc = $nproc;                           // number of files to search into"   >> ${Cpp_code_Folder}/user_inputs_point_2D.h
echo "    float dt = $dt;                                     // timestep"                         >> ${Cpp_code_Folder}/user_inputs_point_2D.h
echo "    float Bx0 = $Bx0;                                   // base field x"                     >> ${Cpp_code_Folder}/user_inputs_point_2D.h
echo "    char output_filename[] = \"${output_filename}\";    // output file"                      >> ${Cpp_code_Folder}/user_inputs_point_2D.h

# compile and execute the c++ code
Script_Folder=`pwd`
cd ${Cpp_code_Folder}
c++ -o process-virtual-sat-point-2D process-virtual-sat-point-2D.cpp
scp ${Cpp_code_Folder}/process-virtual-sat-point-2D ${Input_Data_Folder}/process-virtual-sat-point-2D
cd ${Input_Data_Folder}
./process-virtual-sat-point-2D
mv ${Input_Data_Folder}/BestVirtualSatelliteInfoFoundPoint.txt ${Output_Data_Folder}/BestVirtualSatelliteInfoFoundPoint.txt
mv ${Input_Data_Folder}/${output_filename} ${Output_Data_Folder}/${output_filename}
rm ${Input_Data_Folder}/process-virtual-sat-point-2D
cd ${Script_Folder}

output_file_list=${Input_Data_Folder}"/output_filename_*.txt"

for file in ${output_file_list}; do
    mv "${file}" "${Output_Data_Folder}/."
done

echo
echo 'All done! Now go and have a beer...'
echo

exit
