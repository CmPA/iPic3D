#!/bin/bash
########################################################################
#
#    Author:  A.E.Vapirev, KU Leuven, Afdeling Plasma-astrofysica
#             2011 Sep 16
#
#    Help:    Compiles and executes a c++ script to compute FFTs of a Virtual Satellite data,
#             extracted beforehand with 'process-virtual-satellite-files-point.sh'.
#             After that it plots Log of the extracted FFTs.
#
#    Change log           :
#
########################################################################

########################## define user inputs ##########################

# c++ code folder
Cpp_code_Folder="/home/alexander/Programs/c/read-virt-sat-files-iPIC/fft-cpp-code-point"

# x,y,z coordinates
xx=24
yy=7.5
zz=5

B0x=0.0097                   # from settings.hdf
dt=0.125                     # from settings.hdf
mass_ratio=-256              # take it from settings.hdf - necessary to compute some reference frequencies
istart_record=00 #50     # (array element) beginning of window which we want to FFT - note that in c++ the first array record is [0]
istop_record=11940 #11940   # (array element) ending of window which we want to FFT
frequency_scale_value="wci"
scale_frequency_yes_no="yes" # if 'yes' then the x-axis and all the reference frequencies are scaled with 1/frequency_scale_value

# input/output folders and filenames
Input_Data_Folder="/home/alexander/Desktop/data/tred46/output-data-point-x${xx}y${yy}z${zz}"
Output_Plot_Folder="."
Gnuplot_Script_name="fft-plot-virt-sat-point.gnuplot"
Outout_Plot_Filename="fft-virt-sat-data-plot-point-x${xx}y${yy}z${zz}.jpg"

# varibles
Variables=( "rho_tot_e" ) #( "Bx" "By" "Bz" "Ex" "Ey" "Ez" "V_tot_e" "V_tot_i" "rho_tot_e" "rho_tot_i" )
# Variables are multiplied by the normalizers - for each variable there must be a normalizer
Variables_Normalizers=( "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" )

# complete variable list in the exact output order - it is taken form the corresponding c++ code
Complete_Var_List=( "Bx" "By" "Bz" "Ex" "Ey" "Ez" "Jx_tot_e" "Jy_tot_e" "Jz_tot_e" "Jx_tot_i" "Jy_tot_i" "Jz_tot_i" "rho_tot_e" "rho_tot_i" "Vx_tot_e" "Vy_tot_e" "Vz_tot_e" "Vx_tot_i" "Vy_tot_i" "Vz_tot_i" "B" "E" "V_tot_e" "V_tot_i" "ax_tot_e" "ay_tot_e" "az_tot_e" "ax_tot_i" "ay_tot_i" "az_tot_i" "a_tot_e" "a_tot_i" )

const_w=0.137 #some measured frequency to compare with the rest - this particular is in [rad/wci]

# gnuplot controls

### IMPORTANT: Labels and tick notation change at the bottom plot ###

# stacked plots dimensions for single panel plot
SinglePlotHeight=600
plotwidth=1000

# Define the terminal and total plot size
set_terminal="jpeg size plotwidth,plotheight"

# In order to leave room for axis and tic labels underneath, a N+1 or N+2 -plot layout is created but only use the top N slots.
Nplus="Nplots+1"

# plot title
multiplot_title="FFT of Virt.Sat. Point Data at x=${xx}, y=${yy}, z=${zz} tred46"
font="Helvetica,20"

# margins between up/botom plots
tmargin=0
bmargin=0

# start/stop values for time axis (X) for left (1) and right (2) plots
xstart=0.001
xstop=200.0

# left and right plot alignement
Single_Plot_SizeW="1.0"                               # left plot width - play with this to adjust plot alignment
Single_Plot_SizeH="1.0/(Nplus+1)"                     # single plot height - play with this to adjust plot alignment

# margins from edge of canvas
left_margin=15                                        # play with this to adjust plot alignment
right_margin=3                                        # play with this to adjust plot alignment

# range
y_range_expander=1.0 # if wanted, expand the cbrange option in gnuplot to make sure lowest and highest values are being plotted and not cut out

# ticks
xtics=""
xlabel="Frequency [radians/wci]"
x_format="10^{%L}" #"%3.2e"

# format of y-axis
y_format="%L" #"%3.2e"
min_y=1
max_y=5000

########################### end user inputs ############################


# compute the FFT using a c++ code - compile and execute the c++ code
Script_Folder=`pwd`
cd ${Cpp_code_Folder}
c++ -o fft-virtsat-point-data fft-virtsat-point-data.cpp -lfftw3 -lm
scp ${Cpp_code_Folder}/fft-virtsat-point-data ${Input_Data_Folder}/fft-virtsat-point-data
cd ${Input_Data_Folder}
./fft-virtsat-point-data $mass_ratio $istart_record $istop_record $B0x $dt
rm ${Input_Data_Folder}/fft-virtsat-point-data
cd ${Script_Folder}

# generate the gnuplot code

# automatically find the number of plots from the number of variables
Nplots=${#Variables[@]}

# remove previous gnuplot script
if [ -e ${Output_Plot_Folder}/${Gnuplot_Script_name} ]; then
   rm ${Output_Plot_Folder}/${Gnuplot_Script_name}
fi

# create output folder if not existing
if [ ! -e ${Output_Plot_Folder} ]; then
   mkdir ${Output_Plot_Folder}
fi

# generate gnuplot script

echo "# Stacked XY plots for Virtual Satellite Point Probe"                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "#"                                                                            >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "# If there is bad data points, simply put/replace with a non-numeric character in front and it will skip it." >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "reset"                                                                        >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "############# USER INPUTS #############"                                      >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "### IMPORTANT: Labels and tick notation change at the botom plot ###"         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "# Number of stacked plots"                                                    >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "Nplots=${Nplots}"                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "SinglePlotHeight=${SinglePlotHeight}"                                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "plotwidth=${plotwidth}"                                                       >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "# Define the total plot size"                                                 >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "plotheight=Nplots*SinglePlotHeight"                                           >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set terminal ${set_terminal}"                                                 >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set output \"${Output_Plot_Folder}/${Outout_Plot_Filename}\""                 >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "# In order to leave room for axis and tic labels underneath, a N+1 or N+2 -plot layout is created but only use the top N slots." >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "Nplus=${Nplus}"                                                               >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set multiplot layout Nplus,1 title \"${multiplot_title}\" font \"${font}\""   >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set tmargin ${tmargin}"                                                       >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set bmargin ${bmargin}"                                                       >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "xstart=${xstart}"                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "xstop=${xstop}"                                                               >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "SinglePlotSizeW=${Single_Plot_SizeW}"                                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "SinglePlotSizeH=${Single_Plot_SizeH}"                                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set lmargin ${left_margin}"                                                   >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set rmargin ${right_margin}"                                                  >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "########### END USER INPUTS ###########"                                      >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "# All plots except the bottom one will not have x-axis notation"              >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set format x \"\""                                                            >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set xlabel \"\""                                                              >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set format y \"${y_format}\""                                                 >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set xrange [${xstart}:${xstop}]"                                              >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
#echo "set size SinglePlotSizeW,SinglePlotSizeH"                                     >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set xtics ${xtics}"                                                           >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set log xy"                                                                   >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "## Start form the top - First Figure"                                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}

# read the five referenve frequencies computed with the FFT c++ code
line=`cat ${Input_Data_Folder}/reference_frequencies_fft_sat_output_point.txt`
linearr=( ${line} )
wci=${linearr[0]};
wpi=${linearr[1]};
wce=${linearr[2]};
wlh=${linearr[3]};
wsm=${linearr[4]};

#wci=`echo "$wci * 6.283 "| bc`
#wpi=`echo "$wpi * 6.283 "| bc`
#wce=`echo "$wce * 6.283 "| bc`
#wlh=`echo "$wlh * 6.283 "| bc`
#wsm=`echo "$wsm * 6.283 "| bc`

frequency_axis_scale=1.0
if [ ${scale_frequency_yes_no} = "yes" ]; then
   if [ "$frequency_scale_value" = "wci" ]; then 
      frequency_axis_scale=`echo "scale=3; 1.0 / $wci" | bc`
   fi
   if [ "$frequency_scale_value" = "wpi" ]; then 
      frequency_axis_scale=`echo "scale=3; 1.0 / $wpi" | bc`
   fi
   if [ "$frequency_scale_value" = "wce" ]; then 
      frequency_axis_scale=`echo "scale=3; 1.0 / $wce" | bc`
   fi
   if [ "$frequency_scale_value" = "wlh" ]; then 
      frequency_axis_scale=`echo "scale=3; 1.0 / $wlh" | bc`
   fi
   if [ "$frequency_scale_value" = "wsm" ]; then 
      frequency_axis_scale=`echo "scale=3; 1.0 / $wsm" | bc`
   fi
   wci=`echo "scale=3; $wci * $frequency_axis_scale" | bc`
   wpi=`echo "scale=3; $wpi * $frequency_axis_scale" | bc`
   wce=`echo "scale=3; $wce * $frequency_axis_scale" | bc`
   wlh=`echo "scale=3; $wlh * $frequency_axis_scale" | bc`
   wsm=`echo "scale=3; $wsm * $frequency_axis_scale" | bc`
fi
#frequency_axis_scale=`echo "scale=3; $frequency_axis_scale * 6.283" | bc`

echo ""                                                                      >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set key left top"                                                      >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set parametric"                                                        >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set trange [${min_y}:${max_y}]"                                        >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "const1=${wci}"                                                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "const2=${wpi}"                                                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "const3=${wce}"                                                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "const4=${wlh}"                                                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "const5=${wsm}"                                                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "constw=${const_w}"                                                     >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                      >> ${Output_Plot_Folder}/${Gnuplot_Script_name}

i=0
for icount in ${Variables[@]}; do

    # determine which row to which variable corresponds for gnuplot's 'using' command
    list_item=0
    var_row=0
    for ivar in ${Complete_Var_List[@]}; do
        if [ "${Complete_Var_List[$list_item]}" = "${Variables[$i]}" ]; then
           var_row=$(( $list_item + 1 + 1 ))    # +1 because list_item starts from 0 and then +1 again because 1st column in the data is the X-asis data
        fi
        list_item=$(( $list_item + 1 ))
    done

    var_row_string="\$""${var_row}"

    # generate the separate plots in gnuplot script
    if [ $(( $i + 1 )) -eq ${Nplots} ]; then
       echo "# Last plot with x-axis notation - ### Adjust notation here ###"       >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo ""                                                                      >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "# Set the x-axis notation for the bottom/last plot"                    >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "unset format"                                                          >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "unset xlabel"                                                          >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set xlabel \"${xlabel}\""                                              >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set format x \"${x_format}\""                                          >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set ylabel \"FFT(${Variables[$i]})\""                                  >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set yrange [${min_y}:${max_y}]"                                        >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "plot \"${Input_Data_Folder}/fft_sat_output_point.txt\" every ::2 using (1.0*\$1):(${Variables_Normalizers[$i]}*$var_row_string) notitle with lines, const1,t title \"wci\", const2,t title \"wpi\", const3,t title \"wce\", const4,t title \"wlh\", const5,t title \"wsm\", constw,t title \"w\"" >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo ""                                                                      >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
    else
       echo "set ylabel \"FFT(${Variables[$i]})\""                                  >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set yrange [${min_y}:${max_y}]"                                        >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
#       echo "plot \"${Input_Data_Folder}/fft_sat_output_point.txt\" every ::2 using ($frequency_axis_scale*\$1):(${Variables_Normalizers[$i]}*$var_row_string) notitle with lines, const1,t notitle, const2,t notitle, const3,t notitle, const4,t notitle, const5,t notitle" >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "plot \"${Input_Data_Folder}/fft_sat_output_point.txt\" every ::2 using ($frequency_axis_scale*\$1):(${Variables_Normalizers[$i]}*$var_row_string) notitle with lines, const1,t title \"wci\", const2,t title \"wpi\", const3,t title \"wce\", const4,t title \"wlh\", const5,t title \"wsm\", constw,t title \"w\"" >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo ""                                                                      >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
    fi
    i=$(( $i + 1 ))
done

echo "unset multiplot"                                                              >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "reset"                                                                        >> ${Output_Plot_Folder}/${Gnuplot_Script_name}

echo
echo "---> Plotting..."
gnuplot ${Output_Plot_Folder}/${Gnuplot_Script_name}

echo
echo 'All done! Now go and have a beer...'
echo

exit
