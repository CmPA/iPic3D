#!/bin/bash
########################################################################
#
#    Author:  A.E.Vapirev, KU Leuven, Afdeling Plasma-astrofysica
#             2011 Oct 20
#
#    Help:    This script utilizes a c++ code to extract a user defined fragment
#             of any quantity from an iPIC3D style VTK and then perform and 2D FFT
#             in space of that slice fragment. Lastly, the script autogenerates
#             a gnuplot script and plots the fragment and its FFT. The VTKs must be written
#             using the hdf-to-vtk-converter.cpp by A.E.Vapirev which is based on
#             the original HDFtoVTK.cpp by S.Markidis.
#             Since the FFt is over space, the axes are scaled as discrete wave number
#             over ion skind depth d_i.
#
#    Change log           :
#
########################################################################

########################## define user inputs ##########################

vtk_filename="/media/Elements/Data_and_Runs/iPIC3D_data/tred46/fields_B_cycle14000.vtk"
# cut_plane can be xy, xz, yz
cut_plane="xz"
# var_type can be "scalar" or "vectorX", "vectorY", "vectorZ", "vector_magnitude" for extracting a vector component or magnitude
var_type="vector_magnitude"
slice_pos=13.5
# begin and end coordinates for the portion of the slice in 2D
coord11=10.0
coord12=30.0
coord21=00.0
coord22=15.0

# scale for the FFT plot axes - take it from settings.hdf - note: wci=B0x
B0x=0.0097           # from settings.hdf
dt=0.125             # from settings.hdf
wci=$B0x
# delta_X and delta_Y are probably box size over ( number of satellites - 1 ) in the respective dimension
delta_1=0.3906       #in units d_i : tred46 we have dx=0.3906 ,dy=0.46875, dz=0.39063
delta_2=0.46875
# FFT points according to the number of probes in each dimension - e.g., tred46 we have grid of probes 96x36x24 - the FFT is N/2 ( if odd then (N-1)/2 )
NumberOfStellites_1=48
NumberOfStellites_2=12
scaleFFTaxis_1="2.0*3.1415927/(${NumberOfStellites_1}*${delta_1})"     # sample wave number in dimension 1
scaleFFTaxis_2="2.0*3.1415927/(${NumberOfStellites_2}*${delta_2})"     # sample wave number in dimension 1

# plot axes labeling and title
slice_label_units=""
FFT2D_label_units="k_n/d_i"
Plot_Title="FFT of 2D slice fragment of B-field cycle 14000"

Input_Data_Folder="."
Output_Plot_Folder="."
Cpp_Code_Folder="."
Gnuplot_Script_name="plot-fft2d-slice-of-vtk.gnuplot"
Outout_Plot_Filename="fft-2D-slice-fragment-VTK-${cut_plane}.eps"

########################### end user inputs ############################

# remove previous inputs.h
if [ -e ${Cpp_Code_Folder}/user_inputs.h ]; then
   rm ${Cpp_Code_Folder}/user_inputs.h
fi

# define slice normal
if [ $cut_plane = "xz" ]; then
   slice_normal="y"
fi
if [ $cut_plane = "xy" ]; then
   slice_normal="z"
fi
if [ $cut_plane = "yz" ]; then
   slice_normal="x"
fi

# create inputs.h for the c++ code

echo "    char vtk_filename[]=\"${vtk_filename}\";"                  >> ${Cpp_Code_Folder}/user_inputs.h
echo "    string slice_normal=\"${slice_normal}\";"                  >> ${Cpp_Code_Folder}/user_inputs.h
echo "    string var_type=\"${var_type}\";"                          >> ${Cpp_Code_Folder}/user_inputs.h
echo "    float slice_pos=${slice_pos};"                             >> ${Cpp_Code_Folder}/user_inputs.h
echo "    float coord11=${coord11};"                                 >> ${Cpp_Code_Folder}/user_inputs.h
echo "    float coord12=${coord12};"                                 >> ${Cpp_Code_Folder}/user_inputs.h
echo "    float coord21=${coord21};"                                 >> ${Cpp_Code_Folder}/user_inputs.h
echo "    float coord22=${coord22};"                                 >> ${Cpp_Code_Folder}/user_inputs.h

# remove previous gnuplot script
if [ -e ${Output_Plot_Folder}/${Gnuplot_Script_name} ]; then
   rm ${Output_Plot_Folder}/${Gnuplot_Script_name}
fi

# create output folder if not existing
if [ ! -e ${Output_Plot_Folder} ]; then
   mkdir ${Output_Plot_Folder}
fi

# compile and execute the c++ code

c++ -o fft2d-slice-of-vtk fft2d-slice-of-vtk.cpp -lfftw3 -lm
./fft2d-slice-of-vtk

# read the four indices computed with the c++ code
line=`cat ${Cpp_Code_Folder}/indices.txt`
linearr=( ${line} )
index11=${linearr[0]}
index12=${linearr[1]}
index21=${linearr[2]}
index22=${linearr[3]}

# generate gnuplot

if [ ${slice_normal} = "x" ]; then
   axis1="Y"
   axis2="Z"
fi
if [ ${slice_normal} = "y" ]; then
   axis1="X"
   axis2="Z"
fi
if [ ${slice_normal} = "z" ]; then
   axis1="X"
   axis2="Y"
fi

echo "reset"                                                                                                        >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "############# USER INPUTS #############"                                                                      >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "# Number of stacked plots"                                                                                    >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "Nplots=2"                                                                                                     >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "SinglePlotHeight=400"                                                                                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "plotwidth=1000"                                                                                               >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "# Define the total plot size"                                                                                 >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "plotheight=Nplots*SinglePlotHeight"                                                                           >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set term postscript eps enhanced"                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set output \"${Output_Plot_Folder}/${Outout_Plot_Filename}\""                                                 >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "# In order to leave room for axis and tic labels underneath, a N+1 or N+2 -plot layout is created but only use the top N slots." >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "Nplus=Nplots+1"                                                                                               >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set multiplot layout Nplus,1 title \"${Plot_Title}\" font \"Helvetica,20\""                                   >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set tmargin 0"                                                                                                >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set bmargin 1"                                                                                                >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "leftmargin=5"                                                                                                 >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "rightmargin=5"                                                                                                >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "coord11=$coord11"                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "coord12=$coord12"                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "coord21=$coord21"                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "coord22=$coord22"                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "index11=$index11"                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "index12=$index12"                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "index21=$index21"                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "index22=$index22"                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "scaleX=($coord12-$coord11)/($index12-$index11)"                                                               >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "scaleY=($coord22-$coord21)/($index22-$index21)"                                                               >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "scaleXfft=($index12-$index11)*${scaleFFTaxis_1}"                                                              >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "scaleYfft=($index22-$index21)*${scaleFFTaxis_2}"                                                              >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "########### END USER INPUTS ###########"                                                                      >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set pm3d map"                                                                                                 >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "# All plots except the bottom one will not have x-axis notation"                                              >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "unset format"                                                                                                 >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "## Start form the top - First Figure"                                                                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set xlabel \"$axis1 ${slice_label_units}\""                                                                   >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set ylabel \"$axis2 ${slice_label_units}\""                                                                   >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set cbrange [:]"                                                                                       >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set palette color"                                                                                            >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set palette model RGB"                                                                                        >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set palette defined"                                                                                          >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set format cb \"%g\""                                                                                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set xrange [$index11*scaleX:$index12*scaleX]"                                                                 >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set yrange [$index21*scaleY:$index22*scaleY]"                                                                 >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "splot \"${Output_Plot_Folder}/output_2D_slice_VTK.txt\" matrix using ($coord11+scaleX*\$1):($coord21+scaleY*\$2):(1.0*\$3) notitle" >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set xlabel \"${FFT2D_label_units}\""                                                                          >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set ylabel \"${FFT2D_label_units}\""                                                                          >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set cbrange [:]"                                                                                       >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set palette color"                                                                                            >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set palette model RGB"                                                                                        >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set palette defined"                                                                                          >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set logscale zcb"                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set format cb \"%L^{10}\""                                                                                    >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set xrange [0*scaleXfft:(index12-index11)*scaleXfft]"                                                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set yrange [0*scaleYfft:((index22-index21)/2+1)*scaleYfft]"                                                   >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "splot \"${Output_Plot_Folder}/output_fft_2D_slice_VTK.txt\" matrix using (scaleXfft*\$1):(scaleYfft*\$2):(1.0*\$3) notitle" >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "unset multiplot"                                                                                              >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "reset"                                                                                                        >> ${Output_Plot_Folder}/${Gnuplot_Script_name}


echo
echo "---> Plotting..."
gnuplot ${Output_Plot_Folder}/${Gnuplot_Script_name}

echo
echo 'All done! Now go and have a beer...'
echo

exit
