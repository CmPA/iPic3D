#!/bin/bash
########################################################################
#
#    Author:  A.E.Vapirev, KU Leuven, Afdeling Plasma-astrofysica
#             2011 Sep 16
#
#    Help:
#
#    Change log           :
#
########################################################################

########################## define user inputs ##########################

# x,y,z coordinates
xx=24
yy=7.5
zz=5

# input/output folders and filenames
Input_Data_Folder_1="/home/alexander/Desktop/data/tred46/output-data-point-x${xx}y${yy}z${zz}"
Input_Data_Folder_2="/home/alexander/Desktop/data/tred46.2/output-data-point-x${xx}y${yy}z${zz}"
Output_Plot_Folder="."
Gnuplot_Script_name="plot-virt-sat-point-stitch.gnuplot"
Outout_Plot_Filename="virt-sat-data-plot-point-x${xx}y${yy}z${zz}.jpg"

# varibles
Variables=( "Bx" "By" "Bz" "Ex" "Ey" "Ez" "rho_tot_e" "rho_tot_i" ) #( "Bx" "By" "Bz" "Ex" "Ey" "Ez" "Jx_tot_e" "Jy_tot_e" "Jz_tot_e" "Jx_tot_i" "Jy_tot_i" "Jz_tot_i" "rho_tot_e" "rho_tot_i" "Vx_tot_e" "Vy_tot_e" "Vz_tot_e" "Vx_tot_i" "Vy_tot_i" "Vz_tot_i" "B" "E" "V_tot_e" "V_tot_i" "ax_tot_e" "ay_tot_e" "az_tot_e" "ax_tot_i" "ay_tot_i" "az_tot_i" "a_tot_e" "a_tot_i" )
# Variables are multiplied by the normalizers - for each variable there must be a normalizer
Variables_Normalizers=( "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "-12.5664" "12.5664" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" )

# custom respective variable ranges if wanted - must change code below in 'plot generation section' to use this instead of [minval:maxval]
Variables_Range=( "[0.0004:0.008]" "[-0.001:0.005]" "[-0.0014:0.002]" "[-1e-5:6e-5]" "[-0.00016:2e-5]" "[-6.5e-5:-1.5e-5]" "[0:0.12]" "[0.007:0.017]" "[-0.14:-0.005]" "[0.005:0.14]" )

# complete variable list in the exact output order - it is taken form the corresponding c++ code
Complete_Var_List=( "Bx" "By" "Bz" "Ex" "Ey" "Ez" "Jx_tot_e" "Jy_tot_e" "Jz_tot_e" "Jx_tot_i" "Jy_tot_i" "Jz_tot_i" "rho_tot_e" "rho_tot_i" "Vx_tot_e" "Vy_tot_e" "Vz_tot_e" "Vx_tot_i" "Vy_tot_i" "Vz_tot_i" "B" "E" "V_tot_e" "V_tot_i" "ax_tot_e" "ay_tot_e" "az_tot_e" "ax_tot_i" "ay_tot_i" "az_tot_i" "a_tot_e" "a_tot_i" )

# gnuplot controls

### IMPORTANT: Labels and tick notation change at the bottom plot ###

# stacked plots dimensions for single panel plot
SinglePlotHeight=200
plotwidth=1000

# Define the terminal and total plot size
set_terminal="jpeg size plotwidth,plotheight"

# In order to leave room for axis and tic labels underneath, a N+1 or N+2 -plot layout is created but only use the top N slots.
Nplus="Nplots+1"

# plot title
multiplot_title="Virt.Sat. Point Data at x=${xx}, y=${yy}, z=${zz} tred46+46.2                       "
font="Helvetica,20"

# margins between up/botom plots
tmargin=0
bmargin=0.0

# start/stop values for time axis (X) for left (1) and right (2) plots
xstart1=0
xstop1=12.125
xstart2=0
xstop2=8

# left and right plot alignement
Grid_Ratio_Scaler=0.830                                                               # stitch parameter - play with this to adjust plot alignment
Single_Plot_SizeW_left="1.0*(xstop1-xstart1)/(xstop2+xstop1)*${Grid_Ratio_Scaler}"    # left plot width - play with this to adjust plot alignment
Single_Plot_SizeW_right="1.0*(xstop2-xstart2)/(xstop2+xstop1)*${Grid_Ratio_Scaler}"   # right plot width - play with this to adjust plot alignment
Single_Plot_SizeH="1.0/(Nplus+1)"                                                     # single plot height - play with this to adjust plot alignment

# margins from edge of canvas
left_margin=15                                        # play with this to adjust plot alignment
right_margin=10                                       # play with this to adjust plot alignment

# range
y_range_expander=1.0 # if wanted, expand the cbrange option in gnuplot to make sure lowest and highest values are being plotted and not cut out

# ticks
xtics_1="2"
xtics_2="(\"12000\" 2000, \"14000\" 4000, \"16000\" 6000)"    # use for custom tcs on th esecond plot if wanted
xlabel="                                       Cyclotron periods"

# format of colorbar
y_format="%g" #"%3.2e"

########################### end user inputs ############################

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
echo "set multiplot layout Nplus,2 title \"${multiplot_title}\" font \"${font}\""   >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set tmargin ${tmargin}"                                                       >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set bmargin ${bmargin}"                                                       >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "xstart1=${xstart1}"                                                           >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "xstop1=${xstop1}"                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "xstart2=${xstart2}"                                                           >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "xstop2=${xstop2}"                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "SinglePlotSizeW_left=${Single_Plot_SizeW_left}"                               >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "SinglePlotSizeW_right=${Single_Plot_SizeW_right}"                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "SinglePlotSizeH=${Single_Plot_SizeH}"                                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "leftmargin=${left_margin}"                                                    >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "rightmargin=${right_margin}"                                                  >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "########### END USER INPUTS ###########"                                      >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "# All plots except the bottom one will not have x-axis notation"              >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set format x \"\""                                                            >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set xlabel \"\""                                                              >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "## Start form the top - First Figure"                                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}

i=0
for icount in ${Variables[@]}; do

    # read the min and max for each variable and use it to scale the colorbar for the stitched plots
    line1=`cat ${Input_Data_Folder_1}/output_filename_${Variables[$i]}_minmax.txt`
    linearr1=( ${line1} )
    minval1=${linearr1[0]}
    maxval1=${linearr1[1]}

    line2=`cat ${Input_Data_Folder_2}/output_filename_${Variables[$i]}_minmax.txt`
    linearr2=( ${line2} )
    minval2=${linearr2[0]}
    maxval2=${linearr2[1]}

    resultmin=`expr $minval1 \< $minval2`
    resultmax=`expr $maxval1 \> $maxval2`

    if [ "$resultmin" -eq "1" ]; then
       minval=$minval1
    else
       minval=$minval2
    fi
    if [ "$resultmax" -eq "1" ]; then
       maxval=$maxval1
    else
       maxval=$maxval2
    fi

    maxval=`expr $maxval*$y_range_expander*${Variables_Normalizers[$i]}`
    minval=`expr $minval*$y_range_expander*${Variables_Normalizers[$i]}`

    # Next line checks if ${Variables_Normalizers[$i]} is negative. If true, then Normalizer_sign=1.
    # Then we reverse min and max values on y-scale.
    Normalizer_sign=`echo "${Variables_Normalizers[$i]} < 0" | bc -l`
    if [ $Normalizer_sign -eq 1 ];then
       temp_minval=${maxval}
       temp_maxval=${minval}
       minval=${temp_minval}
       maxval=${temp_maxval}
    fi

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
       echo ""                                                                      >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set format y \"${y_format}\""                                          >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set ylabel \"${Variables[$i]}\""                                       >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set yrange [${minval}:${maxval}]"                                      >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set ytics nomirror"                                                    >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "unset y2tics"                                                          >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set lmargin leftmargin"                                                >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set rmargin 0"                                                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set xrange [${xstart1}:${xstop1}]"                                     >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set size SinglePlotSizeW_left,SinglePlotSizeH"                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set xtics ${xtics_1}"                                                  >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set border 7"                                                          >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "plot \"${Input_Data_Folder_1}/sat_output_point.txt\" every ::2 using (\$1):(${Variables_Normalizers[$i]}*$var_row_string) notitle with lines" >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "unset xlabel"                                                          >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "unset ylabel"                                                          >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set format y \"\""                                                     >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set format y2 \"\""                                                    >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "unset ytics"                                                           >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set y2tics"                                                            >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set lmargin 0"                                                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set rmargin rightmargin"                                               >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set xrange [${xstop1}+${xstart2}:${xstop1}+${xstop2}]"                 >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set size SinglePlotSizeW_right,SinglePlotSizeH"                        >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set border 13"                                                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "plot \"${Input_Data_Folder_2}/sat_output_point.txt\" every ::2 using (${xstop1}+\$1):(${Variables_Normalizers[$i]}*$var_row_string) notitle with lines" >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo ""                                                                      >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
    else
       echo "set ylabel \"${Variables[$i]}\""                                       >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set yrange [${minval}:${maxval}]"                                      >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set format y \"${y_format}\""                                          >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set ytics nomirror"                                                    >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "unset y2tics"                                                          >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set lmargin leftmargin"                                                >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set rmargin 0"                                                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set xrange [${xstart1}:${xstop1}]"                                     >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set size SinglePlotSizeW_left,SinglePlotSizeH"                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set xtics ${xtics_1}"                                                  >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set border 7"                                                          >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "plot \"${Input_Data_Folder_1}/sat_output_point.txt\" every ::2 using (\$1):(${Variables_Normalizers[$i]}*$var_row_string) notitle with lines" >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "unset ylabel"                                                          >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set format y \"\""                                                     >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set format y2 \"\""                                                    >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "unset ytics"                                                           >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set y2tics"                                                            >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set lmargin 0"                                                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set rmargin rightmargin"                                               >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set xrange [${xstop1}+${xstart2}:${xstop1}+${xstop2}]"                 >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set size SinglePlotSizeW_right,SinglePlotSizeH"                        >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set border 13"                                                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "plot \"${Input_Data_Folder_2}/sat_output_point.txt\" every ::2 using (${xstop1}+\$1):(${Variables_Normalizers[$i]}*$var_row_string) notitle with lines" >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
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
