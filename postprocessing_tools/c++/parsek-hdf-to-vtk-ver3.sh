#!/bin/bash

# Author: Alexander E. Vapirev
#         KU Leuven, Afdeling Plasma-astrofysica
#
# Version 3
#
# In version 3, the code works with HDF5 version 1.8.8 and 1.8.9. Previous HDF5 versions are supported in version 2.
#
# HDF5 modificiation by Maria Elena Innocenti, 24 Aug 2012.
#
# This script generates a C++ code which converts Parsek 2D/3D hdf files into vtk files. 
# The backbone is based on the original Stefano Markidis' ConvHDF.cpp code.
# A respective makefile is also automatically created and the script is compiled and executed.
#
# Here the user is given flexible data output and choice of 2D or 3D.
# The user must know what is written in the HDF files. Trying to extract unexisting quantites may not stop the script but will still produce warnings.
# Hint: To see what quantities are written in the procXXXX.hdf files simply type 'h5dump --header procXXXX.hdf' on any 'proc' file.
#
# Usage: Simply define the user options and run the script.
#
# Log: Created May 2011.
#      * Version 1 released
#      16 June 2011 Added compatibility with HDF5 libs newer than 1.6.5.
#      17 June 2011 Added option to write out or skip the mesh points in the VTK file since the grid is structured.
#      18 June 2011 fixed bugs to work with large number of HDF files.
#      28 June 2011 (begin version 2)
#                   Added option to compute and then write in final VTK files some total quantites for ions and electrons.
#                   Major code rearrangement for faster read/write.
#      28 June 2011 Some minor debugging. 
#      * Version 2 released
#      30 June 2011 Correction for the gird spacing dx,dy,dz. The boundaries of the blocks for each proc are being read twice from the HDF
#                   files in each dimension (X,Y,Z) with exception of the block boundaries lying on the surface of the simulation box.
#                   This results in a slighlty larger box size in the VTK files. A respective correction is applied
#                   to the grid spacings, namely, the original [e.g.] dx = Lx/nxc is replaced with dx = Lx/(nxc+XLEN-1)
#                   where XLEN is the number of procs in X-direction.
#       4 July 2011 Some write-out bug fixing.
#                   Added code generator for file Alloc.h (needed for the hdf-to-vtk convertor) in case it does not exist in the same location as this script.
#       5 july 2011 Added normalization option.
#      19 Aug  2011 Input/Output folder handling safety options added.
#      23 Aug  2011 Create a normalization notation file in the 'Output_VTK_folder' for future reference.
#      29 Aug  2011 Option to compute parralel electric field Epar
#      14 Sep  2012 Fixed a bug in double reading of ghost cells.
#      * Version 3 released
#      24 Sep  2012 Support implemented for HDF5 1.8.8. and 1.8.9.
#
# Help: A fairly complete example lists of quantities and their path in the HDFs (separated by empty spaces).
#       Here we separate them into several lists.
#
#       Energy_Var_List="energy/electric energy/kinetic/species_0 energy/kinetic/species_1 energy/magnetic"
#       Fields_Var_List="fields/Bx fields/By fields/Bz fields/Ex fields/Ey fields/Ez"
#       Moments_Var_List="moments/phi moments/rho moments/species_0/Jx moments/species_0/Jy moments/species_0/Jz moments/species_0/pXX moments/species_0/pXY moments/species_0/pXZ moments/species_0/pYY moments/species_0/pYZ moments/species_0/pZZ moments/species_0/rho moments/species_1/Jx moments/species_1/Jy moments/species_1/Jz moments/species_1/pXX moments/species_1/pXY moments/species_1/pXZ moments/species_1/pYY moments/species_1/pYZ moments/species_1/pZZ moments/species_1/rho"
#       Potentials_Var_List="potentials/phi"
#       Vector_Fields_List=( "fields/B" "fields/E" "moments/species_0/J" )
#
#       Important note: If vectors with (x,y,z) components are to be extracted from the HDF files,
#                       the user should put them in the Vector_Fields_List without specifying the separate (x,y,z) components elsewhere.
#                       For example Vector_Fields_List="fields/B" will take care of Bx,By,Bz altogether.
#                       If the (x,y,z) components of a field are specified in another list, then do not put them in the Vector_Fields_List
#                       because this will result in a double declaration for the arrays, files, etc...
#
###################################################################
#######                     User inputs                     #######
###################################################################

# important! - anything with (x,y,z) components can/should be defined in Vector_Fields_List below
#            - if the (x,y,z) components are defined here (anywhere else), they can not be defined in Vector_Fields_List

# general settings
Dimension=3D                       # 2D or 3D topology/simulation
Start_Cycle=0                      # start cycle to read
End_Cycle=15000                    # end cycle to read
delta_cycle_time=1000              # step in [number of cycles] at which each output cycle is written

# list of quantities to be extracted from the hdf - for empty list simply put list equal to ""
Energy_Var_List=""
Fields_Var_List=""
Moments_Var_List="moments/rho moments/species_0/rho moments/species_1/rho moments/species_2/rho moments/species_3/rho"
Potentials_Var_List=""

# prepare VTK files containing these vector quantities for better plotting purposes - for empty list put equal to ""
# this list will dump out the specified fields as vectors in the VTK files
# if the list is not empty, make sure the (x,y,z) components are NOT defined/doubled beforehand (above) - they will be created automatically
Vector_Fields_List="fields/B fields/E moments/species_0/J moments/species_1/J moments/species_2/J moments/species_3/J"

# prepare extra VTK files containing the total density and currents (defined in Moments_Var_List above) and velocity for ions and electrons (density(salar);currents,Vi,Ve(vectors))
Write_Species_Total_Fields="yes"   # {yes,no} - combined/computed quantities for species0+species2 (electrons) and species1+species3 (ions)
                                   # use it only if the proper moments exist in the HDF files and are specified in Moments_Var_List above
                                   # scalars to be computed are "rho_tot_ions rho_tot_electrons"
                                   # vectors to be computed are "J_tot_ions J_tot_electrons V_tot_ions V_tot_electrons"
                                   # velocities are computed as V(x,y,z)=J(x,y,z)/rho

# compute parallel electric field Epar
Compute_Epar="yes"                 # {yes,no} - compute Epar=(ExBx+EyBy+EzBz)/|B|

# normalization - division by a number
Normalize="yes"                                                  # {yes,no}
List_of_Quantities_to_Normalize="rho_tot_ions rho_tot_electrons" # the quantities here must be written exactly as they are in the lists above
List_of_Normalization_Constants="0.0795775 0.0795775"            # respective normalizing factors for each quantity above

# mesh options
Write_Mesh_Points="no"             # {yes,no} - since the grid is structured, the mesh can be omitted in the VTK file to save disk space

# input/output folders
# WARNING: the output VTK folder is removed if it exists and then recreated again.
Input_HDF_folder="/shared/gianni/tred47"           # folder where HDF files reside
Output_VTK_folder="/shared/gianni/tred47/vtk-data" # folder where VTK files will be dumped

# compiler and library options and library path locations
INCLUDE_HDF5_PATH="${HDF5_INC_DIR}/"
LIBRARY_HDF5_PATH="${HDF5_LIB_DIR}/"
HDF5_LIBRARY="-L${LIBRARY_HDF5_PATH} -lhdf5 -L${LIBRARY_HDF5_PATH} -lhdf5_hl"
CPP="g++"
OPTFLAGS="-O2"

# next line may be necessary for the HDF5 libs to work properly - if not, it does not hurt
LD_LIBRARY_PATH=${LIBRARY_HDF5_PATH}:$LD_LIBRARY_PATH

# HDF library version compatibility
export HDF5_DISABLE_VERSION_CHECK=2 # use this line if HDF5 libraries are newer than 1.6.5 - the check can be also set to 1.
                                    # 1 supresses version check but spits out a warning. 2 supresses warnings too.
# additional stuff

# temporary array dimensions
mapx=100
mapy=100
mapz=100

# execute the binary or not after compilation
Execute_Converter_Binary="yes"      # {yes,no}

# compress the output VTK files if wanted (applied only if the converter binary is executed)
compress_vtk_files="no"             # {yes,no}
compressing_program="gzip"          # choose whatever is there
compressing_options=""              # whatever options the program accepts

###################################################################
####### End of user inputs. No need to edit the code below. #######
###################################################################

# Note for the developer: lists with ( ) are used when a loop over the elements is needed.

# Functions

strip_slash(){
# replaces the slash (/) in a string with underscore (_)
# usage: some_string_variable=`strip_slash ${another_string_variable}`
var_no_slash="$1"
slash="/"
underscore="_"
var_no_slash=${var_no_slash//${slash}/${underscore}}
echo $var_no_slash
}

# check if HDF file folder and output VTK folder are the same. If so, force exit - never a good idea to postprocess in the input directory.
if [ ${Input_HDF_folder} = ${Output_VTK_folder} ]; then
   echo
   echo "---> 'Input_HDF_folder' is the same as 'Output_VTK_folder' in 'user inputs'. It is never a good idea to postprocess in the input directory. Fix the path names and rerun the script. Forcing exit..."
   echo
   exit
fi

# create output folder for vtk files if it does not exist
if [ -e ${Output_VTK_folder} ]; then
   echo
   echo "---> The 'Output_VTK_folder' exists. Do you want to remove it and create an empty one? - then press 'y' and 'enter'. If you dont want to delete it, press 'n' and 'enter' and then change it with a new one in 'user inputs'."
   echo
   echo "y/n? and then 'enter'"
   echo
   read yes_no
   if [ "${yes_no}" = "y" ]; then
      rm -r ${Output_VTK_folder}
   else
      echo
      echo "---> exiting..."
      exit
   fi
   echo
fi
mkdir ${Output_VTK_folder}

# create a notation file in the 'Output_VTK_folder' for future reference
if [ "${Normalize}" = "yes" ]; then
   echo "List of the quantities in the VTK files which have been normalized:" >> ${Output_VTK_folder}/NormaliztionExplanation.txt
   echo "${List_of_Quantities_to_Normalize}"                                  >> ${Output_VTK_folder}/NormaliztionExplanation.txt
   echo "List of the respective normalization factors (divide by a factor):"  >> ${Output_VTK_folder}/NormaliztionExplanation.txt
   echo "${List_of_Normalization_Constants}"                                  >> ${Output_VTK_folder}/NormaliztionExplanation.txt
else
   echo "None of the quantities in the VTK files has been normalized."        >> ${Output_VTK_folder}/NormaliztionExplanation.txt
fi

# build the full list of extracted quatities
echo
echo "---> Full list of scalar quantities and vector components to be extracted/computed from HDF to VTK:"
echo
# build a list containing (x,y,z) components of vectors in Vector_Fields_List to be added to Var_List_read to read from HDF files
temp_vector_list=( ${Vector_Fields_List} )
temp_vector_component_list=""
i=0
for ivector in ${temp_vector_list[@]}; do
    temp_vector_component_list=${temp_vector_component_list}" ${temp_vector_list[$i]}x ${temp_vector_list[$i]}y ${temp_vector_list[$i]}z"
    i=$(( $i + 1 ))
done

List_String_Scalar="${Energy_Var_List} ${Fields_Var_List} ${Moments_Var_List} ${Potentials_Var_List}" # all scalars in one list
List_String_Scalar_temp_arrays="${List_String_Scalar}"                                                # used for creation of temporary arrays in Var_list
Var_List_Read=( ${List_String_Scalar} ${temp_vector_component_list} )                                 # used to read the data from HDF
Var_List_All_Scalar_Vector_component="${List_String_Scalar} ${temp_vector_component_list}"            # used to announce on the screen in the beginning what quantities are to be read/written/computed (all of them)
List_VTK_Scalar_String_write=${List_String_Scalar}
List_VTK_Scalar_write=( ${List_VTK_Scalar_String_write} )                                             # used only for writing in vtk files
if [ ${Write_Species_Total_Fields} = "yes" ]; then
   Total_Fields_List_Scalar="rho_tot_ions rho_tot_electrons"
   List_String_Scalar="${List_String_Scalar} ${Total_Fields_List_Scalar}"                             # extend the list with the total densities for ions and electrons
fi
if [ ${Compute_Epar} = "yes" ]; then
   List_String_Scalar="${List_String_Scalar} Epar"                                                    # extend the list with Epar
fi
List_VTK_Scalar_String=${List_String_Scalar}
List_VTK_Scalar=( ${List_VTK_Scalar_String} )                                                         # used for creating vtk files
Var_List=( ${List_String_Scalar_temp_arrays} ${temp_vector_component_list} )                          # used for temporary array definition for all quantites except the computed total ones (if specified)

if [ ${Write_Species_Total_Fields} = "yes" ]; then
   Var_List_All_Scalar_Vector_component="${Var_List_All_Scalar_Vector_component} ${Total_Fields_List_Scalar}"
fi
if [ ${Compute_Epar} = "yes" ]; then
   Var_List_All_Scalar_Vector_component="${Var_List_All_Scalar_Vector_component} Epar"
fi
Var_List_All_Scalar_Vector_Announce=( ${Var_List_All_Scalar_Vector_component} )
i=0
for icount in ${Var_List_All_Scalar_Vector_Announce[@]}; do
    echo "---> $(( $i + 1 )) : ${Var_List_All_Scalar_Vector_Announce[$i]}"
    i=$(( $i + 1 ))
done

echo
echo "---> Full list of vector quantities to be extracted/computed from HDF to VTK:"
echo
List_String_Vector=${Vector_Fields_List}
List_VTK_Vector=( ${List_String_Vector} )                                                             # used only for writing in VTK files at the end
if [ ${Write_Species_Total_Fields} = "yes" ]; then
   Total_Fields_List_Vector="J_tot_ions V_tot_ions J_tot_electrons V_tot_electrons"
   List_String_Vector="${List_String_Vector} ${Total_Fields_List_Vector}"
fi

Vector_List=( ${List_String_Vector} )
i=0
for icount in ${Vector_List[@]}; do
    echo "---> $(( $i + 1 )) : ${Vector_List[$i]}"
    i=$(( $i + 1 ))
done
echo

# normalization lists
Quantities_to_Normalize=( ${List_of_Quantities_to_Normalize} )
Normalization_Constants=( ${List_of_Normalization_Constants} )

# begin coversion from HDF to VTK
echo
if [ ${Dimension} = "2D" ];then
   echo ">>> Begin conversion Parsek2D HDF to VTK <<<"
fi
if [ ${Dimension} = "3D" ];then
   echo ">>> Begin conversion iPIC3D HDF to VTK <<<"
fi
echo

# filename of the C++ hdf-to-vtk converter script
ConverterFileNameBinary="ConvParsek_HDFtoVTK"
ConverterFileName="${ConverterFileNameBinary}.cpp"

# we remove the file ${ConverterFileName} - we need to create a new one
if [ -e ${ConverterFileName} ];then
   rm ${ConverterFileName}
fi
if [ -e ${ConverterFileNameBinary} ]; then
   rm ${ConverterFileNameBinary}
fi

# clean up the old makefile if it exists and create a new one

INC_HDF5="-I${INCLUDE_HDF5_PATH}"
LIB_HDF5="-L${LIBRARY_HDF5_PATH}"
HDF5LIBS="${HDF5_LIBRARY}"

if [ -e makefile.conv ]; then
   rm makefile.conv
fi

echo "#makefile for PARSEK ${Dimension} PROJECT                                                                       " >> makefile.conv
echo "                                                                                                                " >> makefile.conv
echo "${ConverterFileNameBinary}: ${ConverterFileName}                                                                " >> makefile.conv
echo "	${CPP} ${OPTFLAGS} -o ${ConverterFileNameBinary} ${INC_HDF5} ${ConverterFileName} ${LIB_HDF5} ${HDF5LIBS} " >> makefile.conv
echo "                                                                                                                " >> makefile.conv
echo "clean:                                                                                                          " >> makefile.conv
echo "	rm -rf ${ConverterFileNameBinary}                                                                         " >> makefile.conv
echo "                                                                                                                " >> makefile.conv

# Building the main c++ HDF-to-VTK reader/writer

# create the initial part of ConverterFileName which reads settings.hdf

echo "/***************************************************************************                        " >> ${ConverterFileName}
echo "            convHDF5.cpp  -  Convert program to open Parsek Output                                  " >> ${ConverterFileName}
echo "                             -------------------                                                    " >> ${ConverterFileName}
echo "    begin                : Jun 2008                                                                 " >> ${ConverterFileName}
echo "    copyright            : (C) 2004 by Stefano Markidis, Giovanni Lapenta                           " >> ${ConverterFileName}
echo "                                                                                                    " >> ${ConverterFileName}
echo "    Change log           : 2011 May - added more output flexibility via bash - A.E.Vapirev          " >> ${ConverterFileName}
echo "                                                                                                    " >> ${ConverterFileName}
echo " ************************************************************************** */                      " >> ${ConverterFileName}
echo "                                                                                                    " >> ${ConverterFileName}
echo "#include \"hdf5.h\"                                                                                 " >> ${ConverterFileName}
echo "#include \"Alloc.h\"                                                                                " >> ${ConverterFileName}
echo "#include \"math.h\"                                                                                 " >> ${ConverterFileName}
echo "                                                                                                    " >> ${ConverterFileName}
echo "#include <iostream>                                                                                 " >> ${ConverterFileName}
echo "#include <fstream>                                                                                  " >> ${ConverterFileName}
echo "#include <string>                                                                                   " >> ${ConverterFileName}
echo "#include <sstream>                                                                                  " >> ${ConverterFileName}
echo "                                                                                                    " >> ${ConverterFileName}
echo "using std::string;                                                                                  " >> ${ConverterFileName}
echo "using std::stringstream;                                                                            " >> ${ConverterFileName}
echo "using std::ofstream;                                                                                " >> ${ConverterFileName}
echo "using std::cout;                                                                                    " >> ${ConverterFileName}
echo "using std::endl;                                                                                    " >> ${ConverterFileName}
echo "                                                                                                    " >> ${ConverterFileName}
echo "                                                                                                    " >> ${ConverterFileName}
echo "                                                                                                    " >> ${ConverterFileName}
echo "int main (int argc, char **argv) {                                                                  " >> ${ConverterFileName}
echo "  // cycle we want to open                                                                          " >> ${ConverterFileName}
echo "  int n_cycle;                                                                                      " >> ${ConverterFileName}
echo "  sscanf(argv[1],\"%d\",&n_cycle);                                                                  " >> ${ConverterFileName}
echo "  // hdf stuff                                                                                      " >> ${ConverterFileName}
echo "  hid_t    file_id;                                                                                 " >> ${ConverterFileName}
echo "  hid_t    dataset_id;                                                                              " >> ${ConverterFileName}
echo "  herr_t   status;                                                                                  " >> ${ConverterFileName}
echo "  // Open the settings file                                                                         " >> ${ConverterFileName}
echo "  cout << \"  \" << endl;                                                                           " >> ${ConverterFileName}
echo "  cout << \"Open the settings file...  \" << endl;                                                  " >> ${ConverterFileName}
echo "  cout << \"OK\" << endl;                                                                           " >> ${ConverterFileName}
echo "  cout << \"  \" << endl;                                                                           " >> ${ConverterFileName}
echo "  file_id = H5Fopen(\"${Input_HDF_folder}/settings.hdf\", H5F_ACC_RDWR, H5P_DEFAULT);               " >> ${ConverterFileName}
echo "  if (file_id < 0){                                                                                 " >> ${ConverterFileName}
echo "     cout << \"couldn't open file: ${Input_HDF_folder}/settings.hdf\" << endl;                      " >> ${ConverterFileName}
echo "     return -1;                                                                                     " >> ${ConverterFileName}
echo "  }                                                                                                 " >> ${ConverterFileName}
echo "  // First read the topology                                                                        " >> ${ConverterFileName}
echo "  int nproc;                                                                                        " >> ${ConverterFileName}
echo "  dataset_id = H5Dopen(file_id, \"/topology/Nprocs\", H5P_DEFAULT);                                               " >> ${ConverterFileName}
echo "  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&nproc);               " >> ${ConverterFileName}
echo "  status = H5Dclose(dataset_id);                                                                    " >> ${ConverterFileName}
echo "  int XLEN;                                                                                         " >> ${ConverterFileName}
echo "  dataset_id = H5Dopen(file_id, \"/topology/XLEN\", H5P_DEFAULT);                                                " >> ${ConverterFileName}
echo "  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&XLEN);                " >> ${ConverterFileName}
echo "  status = H5Dclose(dataset_id);                                                                    " >> ${ConverterFileName}
echo "  int YLEN;                                                                                         " >> ${ConverterFileName}
echo "  dataset_id = H5Dopen(file_id, \"/topology/YLEN\", H5P_DEFAULT);                                                 " >> ${ConverterFileName}
echo "  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&YLEN);                " >> ${ConverterFileName}
echo "  status = H5Dclose(dataset_id);                                                                    " >> ${ConverterFileName}
if [ ${Dimension} = "3D" ];then
   echo "  int ZLEN;                                                                                      " >> ${ConverterFileName}
   echo "  dataset_id = H5Dopen(file_id, \"/topology/ZLEN\", H5P_DEFAULT);                                              " >> ${ConverterFileName}
   echo "  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&ZLEN);             " >> ${ConverterFileName}
   echo "  status = H5Dclose(dataset_id);                                                                 " >> ${ConverterFileName}
fi
echo "                                                                                                    " >> ${ConverterFileName}
echo "  // read Lx                                                                                        " >> ${ConverterFileName}
echo "  double Lx;                                                                                        " >> ${ConverterFileName}
echo "  dataset_id = H5Dopen(file_id, \"/collective/Lx\", H5P_DEFAULT);                                                " >> ${ConverterFileName}
echo "  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Lx);               " >> ${ConverterFileName}
echo "  status = H5Dclose(dataset_id);                                                                    " >> ${ConverterFileName}
echo "  // read Ly                                                                                        " >> ${ConverterFileName}
echo "  double Ly;	                                                                                  " >> ${ConverterFileName}
echo "  dataset_id = H5Dopen(file_id, \"/collective/Ly\", H5P_DEFAULT);                                                " >> ${ConverterFileName}
echo "  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Ly);               " >> ${ConverterFileName}
echo "  status = H5Dclose(dataset_id);                                                                    " >> ${ConverterFileName}
if [ ${Dimension} = "3D" ];then
   echo "  // read Lz                                                                                     " >> ${ConverterFileName}
   echo "  double Lz;                                                                                     " >> ${ConverterFileName}
   echo "  dataset_id = H5Dopen(file_id, \"/collective/Lz\", H5P_DEFAULT);                                             " >> ${ConverterFileName}
   echo "  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Lz);            " >> ${ConverterFileName}
   echo "  status = H5Dclose(dataset_id);                                                                 " >> ${ConverterFileName}
fi
echo "  // read nxc                                                                                       " >> ${ConverterFileName}
echo "  int nxc;                                                                                          " >> ${ConverterFileName}
echo "  dataset_id = H5Dopen(file_id, \"/collective/Nxc\", H5P_DEFAULT);                                               " >> ${ConverterFileName}
echo "  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&nxc);                 " >> ${ConverterFileName}
echo "  status = H5Dclose(dataset_id);                                                                    " >> ${ConverterFileName}
echo "  // read nyc                                                                                       " >> ${ConverterFileName}
echo "  int nyc; 	                                                                                        " >> ${ConverterFileName}
echo "  dataset_id = H5Dopen(file_id, \"/collective/Nyc\", H5P_DEFAULT);                                               " >> ${ConverterFileName}
echo "  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&nyc);                 " >> ${ConverterFileName}
echo "  status = H5Dclose(dataset_id);                                                                    " >> ${ConverterFileName}
if [ ${Dimension} = "3D" ];then
   echo "  // read nzc                                                                                    " >> ${ConverterFileName}
   echo "  int nzc; 	                                                                                  " >> ${ConverterFileName}
   echo "  dataset_id = H5Dopen(file_id, \"/collective/Nzc\", H5P_DEFAULT);                                            " >> ${ConverterFileName}
   echo "  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&nzc);              " >> ${ConverterFileName}
   echo "  status = H5Dclose(dataset_id);                                                                 " >> ${ConverterFileName}
fi
echo "  // read ns                                                                                        " >> ${ConverterFileName}
echo "  int ns;                                                                                           " >> ${ConverterFileName}
echo "  dataset_id = H5Dopen(file_id, \"/collective/Ns\", H5P_DEFAULT);                                                " >> ${ConverterFileName}
echo "  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&ns);                  " >> ${ConverterFileName}
echo "  // at this point you can close settings                                                           " >> ${ConverterFileName}
echo "  status = H5Fclose(file_id);                                                                       " >> ${ConverterFileName}
echo "  // prepare to read the proc files                                                                 " >> ${ConverterFileName}
echo "  hid_t    *proc_file_id = new hid_t[nproc];                                                        " >> ${ConverterFileName}
echo "  string temp;                                                                                      " >> ${ConverterFileName}
if [ ${Dimension} = "2D" ];then
   echo "  int* cartesian_cor= new int[2];                                                                " >> ${ConverterFileName}
   echo "  int mappa[${mapx}][${mapy}];                                                                   " >> ${ConverterFileName}
fi
if [ ${Dimension} = "3D" ];then
   echo "  int* cartesian_cor= new int[3];                                                                " >> ${ConverterFileName}
   echo "  int mappa[${mapx}][${mapy}][${mapz}];                                                          " >> ${ConverterFileName}
fi
echo "  for (int i=0; i < nproc; i++){                                                                    " >> ${ConverterFileName}
echo "      stringstream ss;                                                                              " >> ${ConverterFileName}
echo "      ss << i;                                                                                      " >> ${ConverterFileName}
echo "      temp = \"${Input_HDF_folder}/proc\" + ss.str() + \".hdf\";                                    " >> ${ConverterFileName}
echo "      proc_file_id[i] = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);                           " >> ${ConverterFileName}
echo "      if (proc_file_id[i] < 0){                                                                     " >> ${ConverterFileName}
echo "         cout << \"couldn't open file:  \"<< temp << endl;                                          " >> ${ConverterFileName}
echo "         return -1;                                                                                 " >> ${ConverterFileName}
echo "      }                                                                                             " >> ${ConverterFileName}
echo "      // read the position in the topology                                                          " >> ${ConverterFileName}
echo "      dataset_id = H5Dopen(proc_file_id[i], \"/topology/cartesian_coord\", H5P_DEFAULT);                         " >> ${ConverterFileName}
echo "      status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,cartesian_cor);    " >> ${ConverterFileName}
if [ ${Dimension} = "2D" ];then
   echo "      mappa[cartesian_cor[0]][cartesian_cor[1]] = i;                                                                                                     " >> ${ConverterFileName}
   echo "      cout << \"file\" << i << \" in topology[\" << cartesian_cor[0] << \"][\" << cartesian_cor[1] << \"]\" << endl;                                     " >> ${ConverterFileName}
fi
if [ ${Dimension} = "3D" ];then
   echo "      mappa[cartesian_cor[0]][cartesian_cor[1]][cartesian_cor[2]] = i;                                                                                   " >> ${ConverterFileName}
   echo "      cout << \"file\" << i << \" in topology[\" << cartesian_cor[0] << \"][\" << cartesian_cor[1] << \"][\" << cartesian_cor[2]  <<\"]\" << endl;       " >> ${ConverterFileName}
fi
echo "      status = H5Dclose(dataset_id);                                                 " >> ${ConverterFileName}
echo "      H5Fclose(proc_file_id[i]);                                                     " >> ${ConverterFileName}
echo "	                                                                               " >> ${ConverterFileName}
echo "  }                                                                                  " >> ${ConverterFileName}
echo "	                                                                               " >> ${ConverterFileName}
echo "  // open the output file                                                            " >> ${ConverterFileName}
echo "  stringstream cc;                                                                   " >> ${ConverterFileName}
echo "  cc << n_cycle;                                                                     " >> ${ConverterFileName}
echo "  // prepare the file                                                                " >> ${ConverterFileName}
if [ ${Dimension} = "2D" ];then
   echo "  int nxn = nxc/XLEN + 1;                                                         " >> ${ConverterFileName}
   echo "  int nyn = nyc/YLEN + 1;                                                         " >> ${ConverterFileName}
   echo "  int nzn = 1;                                                                    " >> ${ConverterFileName}
   echo "  double dx = Lx/(nxc+XLEN-1);                                                    " >> ${ConverterFileName}
   echo "  double dy = Ly/(nyc+YLEN-1);                                                    " >> ${ConverterFileName}
   echo "  double dz = 1;                                                                  " >> ${ConverterFileName}
fi
if [ ${Dimension} = "3D" ];then
   echo "  int nxn = nxc/XLEN;                                                         " >> ${ConverterFileName}
   echo "  int nyn = nyc/YLEN;                                                         " >> ${ConverterFileName}
   echo "  int nzn = nzc/ZLEN;                                                         " >> ${ConverterFileName}
   echo "  double dx = Lx/(nxc);                                                    " >> ${ConverterFileName}
   echo "  double dy = Ly/(nyc);                                                    " >> ${ConverterFileName}
   echo "  double dz = Lz/(nzc);                                                    " >> ${ConverterFileName}
fi
echo "                                                                                     " >> ${ConverterFileName}
echo "  int proc;                                                                          " >> ${ConverterFileName}
echo "  int node;                                                                          " >> ${ConverterFileName}
echo "                                                                                     " >> ${ConverterFileName}

# temporary arrays for vectors
echo "// create temporary storage arrays                                                   " >> ${ConverterFileName}
echo "  double* temp_storageX = new double[(nxn+1)*(nyn+1)*(nzn+1)];                                   " >> ${ConverterFileName}
echo "  double* temp_storageY = new double[(nxn+1)*(nyn+1)*(nzn+1)];                                   " >> ${ConverterFileName}
echo "  double* temp_storageZ = new double[(nxn+1)*(nyn+1)*(nzn+1)];                                   " >> ${ConverterFileName}
if [ ${Dimension} = "2D" ];then
   echo "  double*** vectorX = newArr3(double,(nxn+1)*XLEN,(nyn+1)*YLEN,(nzn+1));                      " >> ${ConverterFileName}
   echo "  double*** vectorY = newArr3(double,(nxn+1)*XLEN,(nyn+1)*YLEN,(nzn+1));                      " >> ${ConverterFileName}
   echo "  double*** vectorZ = newArr3(double,(nxn+1)*XLEN,(nyn+1)*YLEN,(nzn+1));                      " >> ${ConverterFileName}
fi
if [ ${Dimension} = "3D" ];then
   echo "  double*** vectorX = newArr3(double,(nxn+1)*XLEN,(nyn+1)*YLEN,(nzn+1)*ZLEN);                 " >> ${ConverterFileName}
   echo "  double*** vectorY = newArr3(double,(nxn+1)*XLEN,(nyn+1)*YLEN,(nzn+1)*ZLEN);                 " >> ${ConverterFileName}
   echo "  double*** vectorZ = newArr3(double,(nxn+1)*XLEN,(nyn+1)*YLEN,(nzn+1)*ZLEN);                 " >> ${ConverterFileName}
fi
echo "                                                                                     " >> ${ConverterFileName}

# temporary arrays for all scalar variables (and single vector components)
i=0
for icount in ${Var_List[@]}; do
    var_name_underscore=`strip_slash ${Var_List[$i]}`
    echo "  double* temp_storage_${var_name_underscore} = new double[(nxn+1)*(nyn+1)*(nzn+1)];                   " >> ${ConverterFileName}
    if [ ${Dimension} = "2D" ];then
       echo "  double*** vector_${var_name_underscore} = newArr3(double,(nxn+1)*XLEN,(nyn+1)*YLEN,(nzn+1));      " >> ${ConverterFileName}
    fi
    if [ ${Dimension} = "3D" ];then
       echo "  double*** vector_${var_name_underscore} = newArr3(double,(nxn+1)*XLEN,(nyn+1)*YLEN,(nzn+1)*ZLEN); " >> ${ConverterFileName}
    fi
    i=$(( $i + 1 ))
done
echo "                                                                                               " >> ${ConverterFileName}

# temporary variables for the total densities, currents, and velocities
if [ ${Write_Species_Total_Fields} = "yes" ]; then
   # ions
   echo "  double rho_i;          " >> ${ConverterFileName}
   echo "  double j_i_x;          " >> ${ConverterFileName}
   echo "  double j_i_y;          " >> ${ConverterFileName}
   echo "  double j_i_z;          " >> ${ConverterFileName}
   echo "  double v_i_x;          " >> ${ConverterFileName}
   echo "  double v_i_y;          " >> ${ConverterFileName}
   echo "  double v_i_z;          " >> ${ConverterFileName}
   # electrons
   echo "  double rho_e;          " >> ${ConverterFileName}
   echo "  double j_e_x;          " >> ${ConverterFileName}
   echo "  double j_e_y;          " >> ${ConverterFileName}
   echo "  double j_e_z;          " >> ${ConverterFileName}
   echo "  double v_e_x;          " >> ${ConverterFileName}
   echo "  double v_e_y;          " >> ${ConverterFileName}
   echo "  double v_e_z;          " >> ${ConverterFileName}
fi

# temporary variable for Epar
if [ ${Compute_Epar} = "yes" ]; then
   echo "  double E_par;          " >> ${ConverterFileName}
   echo "  double E_dot_B;        " >> ${ConverterFileName}
   echo "  double moduleB;        " >> ${ConverterFileName}
fi

# prepare VTK files and fill in the headers

# initialize VTK files for each scalar (single vector component) variable and in List_VTK_Scalar_String

i=0
for icount in ${List_VTK_Scalar[@]}; do
    var_name_underscore=`strip_slash ${List_VTK_Scalar[$i]}`

    echo "// ${var_name_underscore}                                                    " >> ${ConverterFileName}
    echo "  temp = \"${Output_VTK_folder}/${var_name_underscore}_cycle\"+ cc.str();    " >> ${ConverterFileName}
    echo "  temp += \".vtk\";                                                          " >> ${ConverterFileName}
    echo "  cout << \"Preparing file: \" << temp << endl;                              " >> ${ConverterFileName}
    echo "  ofstream my_file_${var_name_underscore}(temp.c_str());                     " >> ${ConverterFileName}
    # VTK header
    echo "  my_file_${var_name_underscore} << \"# vtk DataFile Version 1.0\" << endl;                           " >> ${ConverterFileName}
    echo "  my_file_${var_name_underscore} << \"${var_name_underscore} Field from Parsek\" << endl;             " >> ${ConverterFileName}
    echo "  my_file_${var_name_underscore} << \"ASCII\" << endl;                                                " >> ${ConverterFileName}
    if [ ${Write_Mesh_Points} = "yes" ]; then
       echo "  my_file_${var_name_underscore} << \"DATASET STRUCTURED_GRID\" << endl;                           " >> ${ConverterFileName}
    else
       echo "  my_file_${var_name_underscore} << \"DATASET STRUCTURED_POINTS\" << endl;                         " >> ${ConverterFileName}
    fi
    if [ ${Dimension} = "2D" ];then
       echo "  my_file_${var_name_underscore} << \"DIMENSIONS \" << nxn*XLEN << \" \" << nyn*YLEN << \" \" << nzn << endl;       " >> ${ConverterFileName}
    fi
    if [ ${Dimension} = "3D" ];then
       echo "  my_file_${var_name_underscore} << \"DIMENSIONS \" << nxn*XLEN << \" \" << nyn*YLEN << \" \" << nzn*ZLEN << endl;  " >> ${ConverterFileName}
    fi
    if [ ${Write_Mesh_Points} = "yes" ]; then
       echo "  my_file_${var_name_underscore} << \"POINTS \" << nxn*nyn*nzn*nproc << \" float\" << endl;        " >> ${ConverterFileName}
    else
       echo "  my_file_${var_name_underscore} << \"ORIGIN 0 0 0 \"  << endl;                                    " >> ${ConverterFileName}
       echo "  my_file_${var_name_underscore} << \"SPACING \" << dx << \" \" << dy << \" \" << dz << endl;      " >> ${ConverterFileName}
    fi

    i=$(( $i + 1 ))
done

# mesh points for scalars/single vector components if specified
if [ ${Write_Mesh_Points} = "yes" ]; then
   echo "  // writing mesh points in VTK for scalars                                                          " >> ${ConverterFileName}
   if [ ${Dimension} = "2D" ];then
      echo "  for (int kk=0; kk < nzn;kk++)                                                                   " >> ${ConverterFileName}
   fi
   if [ ${Dimension} = "3D" ];then
      echo "  for (int kk=0; kk < nzn*ZLEN;kk++)                                                              " >> ${ConverterFileName}
   fi
   echo "    for (int jj=0; jj < nyn*YLEN;jj++)                                                               " >> ${ConverterFileName}
   echo "      for (int ii=0; ii < nxn*XLEN;ii++){                                                            " >> ${ConverterFileName}

   i=0
   for icount in ${List_VTK_Scalar[@]}; do
       var_name_underscore=`strip_slash ${List_VTK_Scalar[$i]}`
       echo "         my_file_${var_name_underscore} << ii*dx  << \" \" << jj*dy  << \" \" << kk*dz << endl;  " >> ${ConverterFileName}
       i=$(( $i + 1 ))
   done

   echo "      }                                                                                              " >> ${ConverterFileName}
fi

# define data type in VTK files for scalars/single vector components
echo "                                                                                                        " >> ${ConverterFileName}
echo "// define data type for scalars/single vector components in VTK files                                   " >> ${ConverterFileName}
echo "                                                                                                        " >> ${ConverterFileName}

i=0
for icount in ${List_VTK_Scalar[@]}; do
    var_name_underscore=`strip_slash ${List_VTK_Scalar[$i]}`
    echo "  my_file_${var_name_underscore} << endl;                                                           " >> ${ConverterFileName}
    echo "  my_file_${var_name_underscore} << \"POINT_DATA \" << nxn*nyn*nzn*nproc << endl;                   " >> ${ConverterFileName}
    echo "  my_file_${var_name_underscore} << \"SCALARS ${var_name_underscore} float 1\" << endl;             " >> ${ConverterFileName}
    echo "  my_file_${var_name_underscore} << \"LOOKUP_TABLE default\" << endl;                               " >> ${ConverterFileName}
    echo "                                                                                                    " >> ${ConverterFileName}
    i=$(( $i + 1 ))
done

# initialize VTK files for each vector variable in Vector_List

i=0
for icount in ${Vector_List[@]}; do
    var_name_underscore=`strip_slash ${Vector_List[$i]}`

    echo "                                                                                                   " >> ${ConverterFileName}
    echo "// Vector ${Vector_List[$i]}                                                                       " >> ${ConverterFileName}
    echo "  temp = \"${Output_VTK_folder}/${var_name_underscore}_cycle\"+ cc.str();                          " >> ${ConverterFileName}
    echo "  temp += \".vtk\";                                                                                " >> ${ConverterFileName}
    echo "  cout << \"Preparing file: \" << temp << endl;                                                    " >> ${ConverterFileName}
    echo "  ofstream my_file_${var_name_underscore}(temp.c_str());                                           " >> ${ConverterFileName}
    echo "  my_file_${var_name_underscore} << \"# vtk DataFile Version 1.0\" << endl;                        " >> ${ConverterFileName}
    # VTK header
    echo "  my_file_${var_name_underscore} << \"${var_name_underscore} Field from Parsek\" << endl;          " >> ${ConverterFileName}
    echo "  my_file_${var_name_underscore} << \"ASCII\" << endl;                                             " >> ${ConverterFileName}
    if [ ${Write_Mesh_Points} = "yes" ]; then
       echo "  my_file_${var_name_underscore} << \"DATASET STRUCTURED_GRID\" << endl;                        " >> ${ConverterFileName}
    else
       echo "  my_file_${var_name_underscore} << \"DATASET STRUCTURED_POINTS\" << endl;                      " >> ${ConverterFileName}
    fi
    if [ ${Dimension} = "2D" ];then
       echo "  my_file_${var_name_underscore} << \"DIMENSIONS \" << nxn*XLEN << \" \" << nyn*YLEN << \" \" << nzn << endl;       " >> ${ConverterFileName}
    fi
    if [ ${Dimension} = "3D" ];then
       echo "  my_file_${var_name_underscore} << \"DIMENSIONS \" << nxn*XLEN << \" \" << nyn*YLEN << \" \" << nzn*ZLEN << endl;  " >> ${ConverterFileName}
    fi
    if [ ${Write_Mesh_Points} = "yes" ]; then
       echo "  my_file_${var_name_underscore} << \"POINTS \" << nxn*nyn*nzn*nproc << \" float\" << endl;        " >> ${ConverterFileName}
    else
       echo "  my_file_${var_name_underscore} << \"ORIGIN 0 0 0 \"  << endl;                                    " >> ${ConverterFileName}
       echo "  my_file_${var_name_underscore} << \"SPACING \" << dx << \" \" << dy << \" \" << dz << endl;      " >> ${ConverterFileName}
    fi

    i=$(( $i + 1 ))
done

# mesh points for vectors if specified
if [ ${Write_Mesh_Points} = "yes" ]; then
   echo "  // writing mesh points in VTK for vectors                                                        " >> ${ConverterFileName}
   if [ ${Dimension} = "2D" ];then
      echo "  for (int kk=0; kk < nzn;kk++)                                                                 " >> ${ConverterFileName}
   fi
   if [ ${Dimension} = "3D" ];then
      echo "  for (int kk=0; kk < nzn*ZLEN;kk++)                                                            " >> ${ConverterFileName}
   fi
   echo "    for (int jj=0; jj < nyn*YLEN;jj++)                                                             " >> ${ConverterFileName}
   echo "      for (int ii=0; ii < nxn*XLEN;ii++){                                                          " >> ${ConverterFileName}

   i=0
   for icount in ${Vector_List[@]}; do
       var_name_underscore=`strip_slash ${Vector_List[$i]}`
       echo "         my_file_${var_name_underscore} << ii*dx  << \" \" << jj*dy  << \" \" << kk*dz << endl;    " >> ${ConverterFileName}
       i=$(( $i + 1 ))
   done

   echo "      }                                                                                                " >> ${ConverterFileName}
fi

# define data type in VTK files for vectors
echo "                                                                                                        " >> ${ConverterFileName}
echo "// define data type for vectors in VTK files                                                            " >> ${ConverterFileName}
echo "                                                                                                        " >> ${ConverterFileName}

i=0
for icount in ${Vector_List[@]}; do
    var_name_underscore=`strip_slash ${Vector_List[$i]}`
    echo "  my_file_${var_name_underscore} << endl;                                                          " >> ${ConverterFileName}
    echo "  my_file_${var_name_underscore} << \"POINT_DATA \" << nxn*nyn*nzn*nproc << endl;                  " >> ${ConverterFileName}
    echo "  my_file_${var_name_underscore} << \"VECTORS ${var_name_underscore} float\" << endl;              " >> ${ConverterFileName}
    echo "                                                                                                   " >> ${ConverterFileName}
    i=$(( $i + 1 ))
done

# reading and writing scalars and single vector components from HDF and writing them in VTK files

# reading the scalars and single vector components

echo "  // Reading all scalars and vector components   " >> ${ConverterFileName}
echo "  cout << \"Reading all scalars and vector components from the HDF dataset\" << endl;           " >> ${ConverterFileName}
echo "  proc = 0;                                      " >> ${ConverterFileName}
echo "  for (int i=0; i < XLEN;i++)                    " >> ${ConverterFileName}
echo "    for (int j=0; j < YLEN;j++)                  " >> ${ConverterFileName}
if [ ${Dimension} = "2D" ];then
    echo "      for (int k=0; k < 1;k++){              " >> ${ConverterFileName}
fi
if [ ${Dimension} = "3D" ];then
    echo "      for (int k=0; k < ZLEN;k++){           " >> ${ConverterFileName}
fi

echo "            stringstream ss;                     " >> ${ConverterFileName}
if [ ${Dimension} = "2D" ];then
   echo "            ss << mappa[i][j];                " >> ${ConverterFileName}
fi
if [ ${Dimension} = "3D" ];then
   echo "            ss << mappa[i][j][k];             " >> ${ConverterFileName}
fi

i=0
for icount in ${Var_List_Read[@]}; do
    var_name_underscore=`strip_slash ${Var_List_Read[$i]}`

    if [ ${Dimension} = "2D" ];then
       echo "                                                                                                " >> ${ConverterFileName}
       echo "            // open file for new data set to read                                               " >> ${ConverterFileName}
       echo "            temp = \"${Input_HDF_folder}/proc\" + ss.str() + \".hdf\";                          " >> ${ConverterFileName}
       echo "            proc_file_id[mappa[i][j]] = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);       " >> ${ConverterFileName}
       echo "            if (proc_file_id[mappa[i][j]] < 0){                                                 " >> ${ConverterFileName}
       echo "               cout << \"couldn't open file:  \"<< temp << endl;                                " >> ${ConverterFileName}
       echo "               return -1;                                                                       " >> ${ConverterFileName}
       echo "            }                                                                                   " >> ${ConverterFileName}
       echo "            temp = \"/${Var_List_Read[$i]}/cycle_\"+ cc.str();                                  " >> ${ConverterFileName}
       echo "            dataset_id = H5Dopen(proc_file_id[mappa[i][j]],temp.c_str(), H5P_DEFAULT);                       " >> ${ConverterFileName}
    fi
    if [ ${Dimension} = "3D" ];then
       echo "                                                                                                " >> ${ConverterFileName}
       echo "            // open file for new data set to read                                               " >> ${ConverterFileName}
       echo "            temp = \"${Input_HDF_folder}/proc\" + ss.str() + \".hdf\";                          " >> ${ConverterFileName}
       echo "            proc_file_id[mappa[i][j][k]] = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);    " >> ${ConverterFileName}
       echo "            if (proc_file_id[mappa[i][j][k]] < 0){                                              " >> ${ConverterFileName}
       echo "               cout << \"couldn't open file:  \"<< temp << endl;                                " >> ${ConverterFileName}
       echo "               return -1;                                                                       " >> ${ConverterFileName}
       echo "            }                                                                                   " >> ${ConverterFileName}
       echo "            temp = \"/${Var_List_Read[$i]}/cycle_\"+ cc.str();                                  " >> ${ConverterFileName}
       echo "            dataset_id = H5Dopen(proc_file_id[mappa[i][j][k]],temp.c_str(), H5P_DEFAULT);                    " >> ${ConverterFileName}
    fi

    echo "            status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storage_${var_name_underscore});          " >> ${ConverterFileName}
    echo "            status = H5Dclose(dataset_id);                                                                                             " >> ${ConverterFileName}
    echo "            node=0;                                                                                                                    " >> ${ConverterFileName}
    echo "            for (int ii=0; ii < (nxn+1);ii++)                                                                                              " >> ${ConverterFileName}
    echo "                for (int jj=0; jj < (nyn+1);jj++)                                                                                          " >> ${ConverterFileName}
    echo "                    for (int kk=0; kk < (nzn+1);kk++){                                                                                     " >> ${ConverterFileName}
    echo "                        vector_${var_name_underscore}[ii + nxn*i][jj + nyn*j][kk + nzn*k] = temp_storage_${var_name_underscore}[node]; " >> ${ConverterFileName}
    echo "                        node++;                                                                                                        " >> ${ConverterFileName}
    echo "                    }                                                                                                                  " >> ${ConverterFileName}
    echo "            // close the file                                                                             " >> ${ConverterFileName}

    if [ ${Dimension} = "2D" ];then
       echo "            H5Fclose(proc_file_id[mappa[i][j]]);                                                       " >> ${ConverterFileName}
    fi
    if [ ${Dimension} = "3D" ];then
       echo "            H5Fclose(proc_file_id[mappa[i][j][k]]);                                                    " >> ${ConverterFileName}
    fi

    i=$(( $i + 1 ))
done

echo "        proc++;            " >> ${ConverterFileName}
echo "      }                    " >> ${ConverterFileName}
echo "                           " >> ${ConverterFileName}

# writing everything in VTK

echo "  // Writing the extracted quantites in VTK                                    " >> ${ConverterFileName}
echo "  cout << \"Writing the extracted quantites in VTK\" << endl;                  " >> ${ConverterFileName}
if [ ${Dimension} = "2D" ];then
   echo "  for (int kk=0; kk < nzn;kk++)                                             " >> ${ConverterFileName}
fi
if [ ${Dimension} = "3D" ];then
   echo "  for (int kk=0; kk < nzn*ZLEN;kk++)                                        " >> ${ConverterFileName}
fi
echo "      for (int jj=0; jj < nyn*YLEN;jj++)                                       " >> ${ConverterFileName}
echo "          for (int ii=0; ii < nxn*XLEN;ii++){                                  " >> ${ConverterFileName}

# writing the scalars and single vector components in VTK
i=0
for icount in ${List_VTK_Scalar_write[@]}; do
    var_name_underscore=`strip_slash ${List_VTK_Scalar_write[$i]}`
    if [ ${Normalize} = "yes" ]; then
        inorm=0
        for quantity in ${Quantities_to_Normalize[@]}; do
            if [ "${Quantities_to_Normalize[$inorm]}" = "${List_VTK_Scalar_write[$i]}" ]; then
               echo "              vector_${var_name_underscore}[ii][jj][kk]=vector_${var_name_underscore}[ii][jj][kk]/${Normalization_Constants[$inorm]};   " >> ${ConverterFileName} 
            fi
            inorm=$(( $inorm + 1 ))
        done
    fi
    echo "              my_file_${var_name_underscore} << vector_${var_name_underscore}[ii][jj][kk] << endl;      " >> ${ConverterFileName}
    i=$(( $i + 1 ))
done

# writing the vectors in VTK
i=0
for icount in ${List_VTK_Vector[@]}; do
    var_name_underscore=`strip_slash ${List_VTK_Vector[$i]}`
    if [ ${Normalize} = "yes" ]; then
        inorm=0
        for quantity in ${Quantities_to_Normalize[@]}; do
            if [ "${Quantities_to_Normalize[$inorm]}" = "${List_VTK_Vector[$i]}" ]; then
               echo "              vector_${var_name_underscore}x[ii][jj][kk]=vector_${var_name_underscore}x[ii][jj][kk]/${Normalization_Constants[$inorm]};   " >> ${ConverterFileName} 
               echo "              vector_${var_name_underscore}y[ii][jj][kk]=vector_${var_name_underscore}y[ii][jj][kk]/${Normalization_Constants[$inorm]};   " >> ${ConverterFileName} 
               echo "              vector_${var_name_underscore}z[ii][jj][kk]=vector_${var_name_underscore}z[ii][jj][kk]/${Normalization_Constants[$inorm]};   " >> ${ConverterFileName} 
            fi
            inorm=$(( $inorm + 1 ))
        done
    fi
    echo "              my_file_${var_name_underscore} << vector_${var_name_underscore}x[ii][jj][kk] << \" \" << vector_${var_name_underscore}y[ii][jj][kk] << \" \" << vector_${var_name_underscore}z[ii][jj][kk] << endl;   " >> ${ConverterFileName} 
    i=$(( $i + 1 ))
done

# writing the total quantities (rho,J,V) for ions and electrons if specified
if [ ${Write_Species_Total_Fields} = "yes" ]; then

   # ions
   echo "              // rho,J,V for ions                                                                             " >> ${ConverterFileName}
   echo "              rho_i=vector_moments_species_1_rho[ii][jj][kk]+vector_moments_species_3_rho[ii][jj][kk];        " >> ${ConverterFileName}
   # normalize ion density
   if [ ${Normalize} = "yes" ]; then
       inorm=0
       for quantity in ${Quantities_to_Normalize[@]}; do
           if [ "${Quantities_to_Normalize[$inorm]}" = "rho_tot_ions" ]; then
              echo "              rho_i=rho_i/${Normalization_Constants[$inorm]};   " >> ${ConverterFileName} 
           fi
           inorm=$(( $inorm + 1 ))
       done
   fi
   echo "              my_file_rho_tot_ions << rho_i << endl;                                                          " >> ${ConverterFileName}

   echo "              j_i_x=vector_moments_species_1_Jx[ii][jj][kk]+vector_moments_species_3_Jx[ii][jj][kk];          " >> ${ConverterFileName}
   echo "              j_i_y=vector_moments_species_1_Jy[ii][jj][kk]+vector_moments_species_3_Jy[ii][jj][kk];          " >> ${ConverterFileName}
   echo "              j_i_z=vector_moments_species_1_Jz[ii][jj][kk]+vector_moments_species_3_Jz[ii][jj][kk];          " >> ${ConverterFileName}
   # normalize ion current
   if [ ${Normalize} = "yes" ]; then
       inorm=0
       for quantity in ${Quantities_to_Normalize[@]}; do
           if [ "${Quantities_to_Normalize[$inorm]}" = "J_tot_ions" ]; then
              echo "              j_i_x=j_i_x/${Normalization_Constants[$inorm]};   " >> ${ConverterFileName} 
              echo "              j_i_y=j_i_y/${Normalization_Constants[$inorm]};   " >> ${ConverterFileName} 
              echo "              j_i_z=j_i_z/${Normalization_Constants[$inorm]};   " >> ${ConverterFileName} 
           fi
           inorm=$(( $inorm + 1 ))
       done
   fi
   echo "              my_file_J_tot_ions << j_i_x << \" \" << j_i_y << \" \" << j_i_z << endl;                        " >> ${ConverterFileName}

   echo "              v_i_x=j_i_x/rho_i;                                                                              " >> ${ConverterFileName}
   echo "              v_i_y=j_i_y/rho_i;                                                                              " >> ${ConverterFileName}
   echo "              v_i_z=j_i_z/rho_i;                                                                              " >> ${ConverterFileName}
   # normalize ion velocity
   if [ ${Normalize} = "yes" ]; then
       inorm=0
       for quantity in ${Quantities_to_Normalize[@]}; do
           if [ "${Quantities_to_Normalize[$inorm]}" = "V_tot_ions" ]; then
              echo "              v_i_x=v_i_x/${Normalization_Constants[$inorm]};   " >> ${ConverterFileName} 
              echo "              v_i_y=v_i_y/${Normalization_Constants[$inorm]};   " >> ${ConverterFileName} 
              echo "              v_i_z=v_i_z/${Normalization_Constants[$inorm]};   " >> ${ConverterFileName} 
           fi
           inorm=$(( $inorm + 1 ))
       done
   fi
   echo "              my_file_V_tot_ions << v_i_x << \" \" << v_i_y << \" \" << v_i_z << endl;                        " >> ${ConverterFileName}

   # electrons
   echo "              // rho,J,V for electrons                                                                        " >> ${ConverterFileName}
   echo "              rho_e=vector_moments_species_0_rho[ii][jj][kk]+vector_moments_species_2_rho[ii][jj][kk];        " >> ${ConverterFileName}
   # normalize electron density
   if [ ${Normalize} = "yes" ]; then
       inorm=0
       for quantity in ${Quantities_to_Normalize[@]}; do
           if [ "${Quantities_to_Normalize[$inorm]}" = "rho_tot_electrons" ]; then
              echo "              rho_e=rho_e/${Normalization_Constants[$inorm]};   " >> ${ConverterFileName} 
           fi
           inorm=$(( $inorm + 1 ))
       done
   fi
   echo "              my_file_rho_tot_electrons << rho_e << endl;                                                     " >> ${ConverterFileName}

   echo "              j_e_x=vector_moments_species_0_Jx[ii][jj][kk]+vector_moments_species_2_Jx[ii][jj][kk];          " >> ${ConverterFileName}
   echo "              j_e_y=vector_moments_species_0_Jy[ii][jj][kk]+vector_moments_species_2_Jy[ii][jj][kk];          " >> ${ConverterFileName}
   echo "              j_e_z=vector_moments_species_0_Jz[ii][jj][kk]+vector_moments_species_2_Jz[ii][jj][kk];          " >> ${ConverterFileName}
   # normalize electron current
   if [ ${Normalize} = "yes" ]; then
       inorm=0
       for quantity in ${Quantities_to_Normalize[@]}; do
           if [ "${Quantities_to_Normalize[$inorm]}" = "J_tot_electrons" ]; then
              echo "              j_e_x=j_e_x/${Normalization_Constants[$inorm]};   " >> ${ConverterFileName} 
              echo "              j_e_y=j_e_y/${Normalization_Constants[$inorm]};   " >> ${ConverterFileName} 
              echo "              j_e_z=j_e_z/${Normalization_Constants[$inorm]};   " >> ${ConverterFileName} 
           fi
           inorm=$(( $inorm + 1 ))
       done
   fi
   echo "              my_file_J_tot_electrons << j_e_x << \" \" << j_e_y << \" \" << j_e_z << endl;                   " >> ${ConverterFileName}

   echo "              v_e_x=j_e_x/rho_e;                                                                              " >> ${ConverterFileName}
   echo "              v_e_y=j_e_y/rho_e;                                                                              " >> ${ConverterFileName}
   echo "              v_e_z=j_e_z/rho_e;                                                                              " >> ${ConverterFileName}
   # normalize electron velocity
   if [ ${Normalize} = "yes" ]; then
       inorm=0
       for quantity in ${Quantities_to_Normalize[@]}; do
           if [ "${Quantities_to_Normalize[$inorm]}" = "V_tot_electrons" ]; then
              echo "              v_e_x=v_e_x/${Normalization_Constants[$inorm]};   " >> ${ConverterFileName} 
              echo "              v_e_y=v_e_y/${Normalization_Constants[$inorm]};   " >> ${ConverterFileName} 
              echo "              v_e_z=v_e_z/${Normalization_Constants[$inorm]};   " >> ${ConverterFileName} 
           fi
           inorm=$(( $inorm + 1 ))
       done
   fi
   echo "              my_file_V_tot_electrons << v_e_x << \" \" << v_e_y << \" \" << v_e_z << endl;                   " >> ${ConverterFileName}
fi

# write Epar in VTK file if specified
if [ ${Compute_Epar} = "yes" ]; then
   echo "              // compute Epar                                                                      " >> ${ConverterFileName}
   echo "              E_dot_B=vector_fields_Ex[ii][jj][kk]*vector_fields_Bx[ii][jj][kk]+vector_fields_Ey[ii][jj][kk]*vector_fields_By[ii][jj][kk]+vector_fields_Ez[ii][jj][kk]*vector_fields_Bz[ii][jj][kk];          " >> ${ConverterFileName}
   echo "              moduleB=sqrt(vector_fields_Bx[ii][jj][kk]*vector_fields_Bx[ii][jj][kk]+vector_fields_By[ii][jj][kk]*vector_fields_By[ii][jj][kk]+vector_fields_Bz[ii][jj][kk]*vector_fields_Bz[ii][jj][kk]);    " >> ${ConverterFileName}
   echo "              E_par=E_dot_B/moduleB;                                                               " >> ${ConverterFileName}
   # normalize Epar
   if [ ${Normalize} = "yes" ]; then
       inorm=0
       for quantity in ${Quantities_to_Normalize[@]}; do
           if [ "${Quantities_to_Normalize[$inorm]}" = "Epar" ]; then
              echo "              E_par=E_par/${Normalization_Constants[$inorm]};   " >> ${ConverterFileName} 
           fi
           inorm=$(( $inorm + 1 ))
       done
   fi
   echo "              my_file_Epar << E_par << endl;                                                       " >> ${ConverterFileName}
fi

echo "          }                                                                       " >> ${ConverterFileName}
echo "                                                                                  " >> ${ConverterFileName}


# finilize the main C++ program

echo "// finilize the main program                           " >> ${ConverterFileName}
echo "                                                       " >> ${ConverterFileName}
echo "  delete[] proc_file_id;                               " >> ${ConverterFileName}
echo "  delete[] temp_storageX;                              " >> ${ConverterFileName}
echo "  delete[] temp_storageY;                              " >> ${ConverterFileName}
echo "  delete[] temp_storageZ;                              " >> ${ConverterFileName}
echo "  delArr3(vectorX,(nxn+1)*XLEN,(nyn+1)*YLEN);                  " >> ${ConverterFileName}
echo "  delArr3(vectorY,(nxn+1)*XLEN,(nyn+1)*YLEN);                  " >> ${ConverterFileName}
echo "  delArr3(vectorZ,(nxn+1)*XLEN,(nyn+1)*YLEN);                  " >> ${ConverterFileName}
echo "                                                       " >> ${ConverterFileName}

i=0
for icount in ${Var_List[@]}; do
    var_name_underscore=`strip_slash ${Var_List[$i]}`
    echo "  delete[] temp_storage_${var_name_underscore};                        " >> ${ConverterFileName}
    # cleaning of the arrays in 2D and 3D is the same because in 2D we still have a (uniform) 3rd dimension
    if [ ${Dimension} = "2D" ];then
       echo "  delArr3(vector_${var_name_underscore},(nxn+1)*XLEN,(nyn+1)*YLEN);         " >> ${ConverterFileName}
    fi
    if [ ${Dimension} = "3D" ];then
       echo "  delArr3(vector_${var_name_underscore},(nxn+1)*XLEN,(nyn+1)*YLEN);         " >> ${ConverterFileName}
    fi
    i=$(( $i + 1 ))
done
echo "                                                                           " >> ${ConverterFileName}

# close the vector VTK files
i=0
for icount in ${Vector_List[@]}; do
    var_name_underscore=`strip_slash ${Vector_List[$i]}`
    echo "// close file for ${var_name_underscore}           " >> ${ConverterFileName}
    echo "  my_file_${var_name_underscore}.close();          " >> ${ConverterFileName}
    i=$(( $i + 1 ))
done

# close the scalar VTK files
i=0
for icount in ${List_VTK_Scalar[@]}; do
    var_name_underscore=`strip_slash ${List_VTK_Scalar[$i]}`
    echo "// close file for ${List_VTK_Scalar[$i]}           " >> ${ConverterFileName}
    echo "  my_file_${var_name_underscore}.close();          " >> ${ConverterFileName}
    i=$(( $i + 1 ))
done
echo "                                           " >> ${ConverterFileName}
echo "  return(0);                               " >> ${ConverterFileName}
echo "}                                          " >> ${ConverterFileName}
echo "// end main program                        " >> ${ConverterFileName}
echo "                                           " >> ${ConverterFileName}


# End of the generated c++ code for the hdf-to-vtk converter

# Begin generate c++ code for file Alloc.h (if it does not exist) which is used by the converter
if [ ! -e Alloc.h ]; then
    echo '  * File Alloc.h does not exist -> will create it.'
    echo '/***************************************************************************'           >> Alloc.h
    echo ''                                                                                       >> Alloc.h
    echo '                             -------------------'                                       >> Alloc.h
    echo '    begin                : Fri Jun 4 2004'                                              >> Alloc.h
    echo '    copyright            : (C) 2004 Los Alamos National Laboratory'                     >> Alloc.h
    echo '    developers           : Stefano Markidis, Giovanni Lapenta'                          >> Alloc.h
    echo '    email                : markidis@lanl.gov, lapenta@lanl.gov'                         >> Alloc.h
    echo ' ***************************************************************************/'          >> Alloc.h
    echo '#ifndef Alloc_H'                                                                        >> Alloc.h
    echo '#define Alloc_H'                                                                        >> Alloc.h
    echo '/** subroutines for allocation and deallocation of arrays 2D, 3D, 4D */'                >> Alloc.h
    echo '/**'                                                                                    >> Alloc.h
    echo ' 2 dimensional arrays'                                                                  >> Alloc.h
    echo '*/'                                                                                     >> Alloc.h
    echo ''                                                                                       >> Alloc.h
    echo '/** The allocator for 2D array */'                                                      >> Alloc.h
    echo 'template <class type>'                                                                  >> Alloc.h
    echo 'inline type **_new_2_array(int sz1,int sz2,type *stupid)'                               >> Alloc.h
    echo '{'                                                                                      >> Alloc.h
    echo ' type **foo;'                                                                           >> Alloc.h
    echo ' //foo = new (type *)[sz1];'                                                            >> Alloc.h
    echo ' foo = new type *[sz1];'                                                                >> Alloc.h
    echo ' for (int i=0;i<sz1;i++) foo[i] = new type[sz2];'                                       >> Alloc.h
    echo ' return foo;'                                                                           >> Alloc.h
    echo '}'                                                                                      >> Alloc.h
    echo '/** macro for allocate 2D array */'                                                     >> Alloc.h
    echo '#define newArr(type,sz1,sz2) _new_2_array((sz1),(sz2),(type *) NULL)'                   >> Alloc.h
    echo ''                                                                                       >> Alloc.h
    echo '/** deallocator for a 2D array*/'                                                       >> Alloc.h
    echo 'template <class type>'                                                                  >> Alloc.h
    echo 'inline void delArr(type **foo,int sz1)'                                                 >> Alloc.h
    echo '{ for (int i=0;i<sz1;i++) delete[] foo[i]; delete[] foo; }'                             >> Alloc.h
    echo ''                                                                                       >> Alloc.h
    echo '/**'                                                                                    >> Alloc.h
    echo ' 3 dimensional arrays'                                                                  >> Alloc.h
    echo '*/'                                                                                     >> Alloc.h
    echo ''                                                                                       >> Alloc.h
    echo '/** The allocator for 3D array */'                                                      >> Alloc.h
    echo 'template <class type>'                                                                  >> Alloc.h
    echo 'inline type ***_new_3_array(int sz1,int sz2,int sz3,type *stupid)'                      >> Alloc.h
    echo '{'                                                                                      >> Alloc.h
    echo ' type ***foo;'                                                                          >> Alloc.h
    echo ' foo = new type **[sz1];'                                                               >> Alloc.h
    echo ' for (int i=0;i<sz1;i++) foo[i] = newArr(type,sz2,sz3);'                                >> Alloc.h
    echo ' return foo;'                                                                           >> Alloc.h
    echo '}'                                                                                      >> Alloc.h
    echo '/** macro for allocate 3D array */'                                                     >> Alloc.h
    echo '#define newArr3(type,sz1,sz2,sz3) _new_3_array((sz1),(sz2),(sz3),(type *) NULL)'        >> Alloc.h
    echo ''                                                                                       >> Alloc.h
    echo '/** deallocator for a 3D array*/'                                                       >> Alloc.h
    echo 'template <class type>'                                                                  >> Alloc.h
    echo 'inline void delArr3(type ***foo,int sz1,int sz2)'                                       >> Alloc.h
    echo '{ for (int i=0;i<sz1;i++) delArr(foo[i],sz2); delete[] foo; }'                          >> Alloc.h
    echo ''                                                                                       >> Alloc.h
    echo '/**'                                                                                    >> Alloc.h
    echo '4 dimensional arrays'                                                                   >> Alloc.h
    echo '*/'                                                                                     >> Alloc.h
    echo ''                                                                                       >> Alloc.h
    echo '/** The allocator for 4D array */'                                                      >> Alloc.h
    echo 'template <class type>'                                                                  >> Alloc.h
    echo 'inline type ****_new_4_array(int sz1,int sz2,int sz3,int sz4,type *stupid)'             >> Alloc.h
    echo '{'                                                                                      >> Alloc.h
    echo ' type ****foo;'                                                                         >> Alloc.h
    echo ' foo = new type ***[sz1];'                                                              >> Alloc.h
    echo ' for (int i=0;i<sz1;i++) foo[i] = newArr3(type,sz2,sz3,sz4);'                           >> Alloc.h
    echo ' return foo;'                                                                           >> Alloc.h
    echo '}'                                                                                      >> Alloc.h
    echo '/** macro for allocate 4D array */'                                                     >> Alloc.h
    echo '#define newArr4(type,sz1,sz2,sz3,sz4) \'                                                >> Alloc.h
    echo '_new_4_array((sz1),(sz2),(sz3),(sz4),(type *) NULL);'                                   >> Alloc.h
    echo ''                                                                                       >> Alloc.h
    echo '/** deallocator for a 4D array*/'                                                       >> Alloc.h
    echo 'template <class type>'                                                                  >> Alloc.h
    echo 'inline void delArr4(type ****foo,int sz1,int sz2,int sz3)'                              >> Alloc.h
    echo '{ for (int i=0;i<sz1;i++) delArr3(foo[i],sz2,sz3); delete[] foo; }'                     >> Alloc.h
    echo ''                                                                                       >> Alloc.h
    echo '#endif'                                                                                 >> Alloc.h
fi
# End of the generated c++ code for file Alloc.h

# compile the converter
echo
echo "  * Compiling..."
echo
make -f makefile.conv
echo

# execute the converter binary if specified

if [ ${Execute_Converter_Binary} = "yes" ]; then
   echo "  * Running the compiled binary..."
   echo
   echo
   # Loop over the cycles - building the time sequence
   num_total_cycles=$(( ( $End_Cycle - $Start_Cycle) / $delta_cycle_time + 1 ))
   for (( icycle=1; icycle<=$num_total_cycles; icycle++ )); do
       cycle_number=$(( $Start_Cycle + ( $icycle - 1 ) * $delta_cycle_time ))
       echo "-------> Processing cycle $cycle_number ..."
       # execute the converter binary
       ./${ConverterFileNameBinary} $cycle_number
   done # end loop over cycle_number
   # compress if specified
   if [ ${compress_vtk_files} = "yes" ]; then
      echo
      Final_Output_VTK_File_List=${Output_VTK_folder}"/*.vtk"
      for i in ${Final_Output_VTK_File_List};do
          current_vtk_file_to_compress=${i}
          echo "   ---> Compressing VTK file ${current_vtk_file_to_compress} ..."
          ${compressing_program} ${compressing_options} ${current_vtk_file_to_compress}
      done
      echo
   fi
else
   echo "---> Execute_Converter_Binary option is set to \"no\". You will have to run it yourself with a proper cycle number, e.g., ./${ConverterFileNameBinary} 1000"
fi

echo
if [ ${Dimension} = "2D" ];then
   echo ">>> End conversion Parsek2D HDF to VTK <<<"
fi
if [ ${Dimension} = "3D" ];then
   echo ">>> End conversion iPIC3D HDF to VTK <<<"
fi
echo
echo 'All done! Now go and have a beer...'
echo

exit

