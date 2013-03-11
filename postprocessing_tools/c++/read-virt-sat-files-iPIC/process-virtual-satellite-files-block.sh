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
XYZ_block="XY"
x_satellite_user=13.0859                # x-coordinate of virtual satellite
y_satellite_user=6.44532                # y-coordinate of virtual satellite
z_satellite_user=0.195311               # z-coordinate of virtual satellite
nproc=3071                              # number of files to search into

# inputs to write out the data for the satellite
Bx0=0.0097                              # base field x
dt=0.125                                # timestep
output_filename="sat_output_block.txt"  # output file

# input and output folder
Cpp_code_Folder="/home/alexander/Programs/c/read-virt-sat-files-iPIC/cpp-code-block"
Input_Data_Folder="/home/alexander/Programs/c/read-virt-sat-files-iPIC/virt-sat-data-folder"
Output_Data_Folder="/home/alexander/Programs/c/read-virt-sat-files-iPIC/output-data-block"

########################### end user inputs ############################

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
if [ -e ${Cpp_code_Folder}/user_inputs_block.h ];then
   rm ${Cpp_code_Folder}/user_inputs_block.h
fi

echo

# create user_inputs_block.h file
echo "    string XYZ_block = \"${XYZ_block}\";                // dimension to scan for satellites"   >> ${Cpp_code_Folder}/user_inputs_block.h
echo "    float x_satellite_user = $x_satellite_user;         // x-coordinate"                       >> ${Cpp_code_Folder}/user_inputs_block.h
echo "    float y_satellite_user = $y_satellite_user;         // y-coordinate"                       >> ${Cpp_code_Folder}/user_inputs_block.h
echo "    float z_satellite_user = $z_satellite_user;         // z-coordinate"                       >> ${Cpp_code_Folder}/user_inputs_block.h
echo "    int const nproc = $nproc;                           // number of files to search into"     >> ${Cpp_code_Folder}/user_inputs_block.h
echo "    float dt = $dt;                                     // timestep"                           >> ${Cpp_code_Folder}/user_inputs_block.h
echo "    float Bx0 = $Bx0;                                   // base field x"                       >> ${Cpp_code_Folder}/user_inputs_block.h
echo "    char output_filename[] = \"${output_filename}\";    // output file"                        >> ${Cpp_code_Folder}/user_inputs_block.h

# compile and execute the c++ code
Script_Folder=`pwd`
cd ${Cpp_code_Folder}
c++ -o process-virtual-sat-block process-virtual-sat-block.cpp
scp ${Cpp_code_Folder}/process-virtual-sat-block ${Input_Data_Folder}/process-virtual-sat-block
cd ${Input_Data_Folder}
./process-virtual-sat-block
mv ${Input_Data_Folder}/BestVirtualSatelliteInfoFoundBlock.txt ${Output_Data_Folder}/BestVirtualSatelliteInfoFoundBlock.txt
mv ${Input_Data_Folder}/${output_filename} ${Output_Data_Folder}/${output_filename}
rm ${Input_Data_Folder}/process-virtual-sat-block
cd ${Script_Folder}

echo
echo 'All done! Now go and have a beer...'
echo

exit
