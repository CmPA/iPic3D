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

# traces
TR_coord=x
TR=13

# input/output folders and filenames
Input_Data_Folder="/home/alexander/data/test/virt-sat-data-folder/temp-trace-x13/output-data-block-single-var-${TR_coord}${TR}.2"
Output_Plot_Folder="."
Gnuplot_Script_name="plot-virt-sat-block-single-var-2D.gnuplot"
Outout_Plot_Filename="virt-sat-data-plot-trace-${TR_coord}${TR}.jpg"

# varibles
Variables=( "Bx" "By" "Bz" "Ex" "Ey" "Ez" "rho_tot_e" "rho_tot_i" )
# Variables are multiplied by the normalizers - for each variable there must be a normalizer
Variables_Normalizers=( "1.0" "1.0" "1.0" "1.0" "1.0" "1.0" "-12.5664" "12.5664" )

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
multiplot_title="Virt.Sat. Trace Data at ${TR_coord}=${TR} run157"
font="Helvetica,20"

# margins between up/botom plots
tmargin=0
bmargin=0

# start/stop values for time axis (X)
xstart=0
xstop=40000

# left and right plot alignement
Single_Plot_SizeW="1.0"                              # left plot width - play with this to adjust plot alignment
Single_Plot_SizeH="1.0/Nplus*1.3"                    # single plot height - play with this to adjust plot alignment

# margins from edge of canvas
left_margin=5                                        # play with this to adjust plot alignment
right_margin=5                                       # play with this to adjust plot alignment

# range
yrange="[1:12]"             # max is number of satellites in the respactive trace direction
ytics="" #"(\"0\" 1, \"2\" 4.8, \"4\" 9.6, \"6\" 14.4, \"8\" 19.2, \"10\" 24)"
ylabel="Z"                  # along this coordinate are the satellite probes
colorbar_range_expander=1.0 # if wanted, expand the cbrange option in gnuplot to make sure lowest and highest values are being plotted and not cut out

# palette
palette_color="color"
palette_model="RGB"
palette="defined"

# ticks
xtics="2000"
xlabel="Cycle"

# format of colorbar
colorbar_format="%g" #"%3.2e"

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

echo "# Stacked pm3d plots for Virtual Satellite Traces"                            >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
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
echo "leftmargin=${left_margin}"                                                    >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "rightmargin=${right_margin}"                                                  >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set yrange ${yrange}"                                                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set palette ${palette_color}"                                                 >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set palette model ${palette_model}"                                           >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set palette ${palette}"                                                       >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set format cb \"${colorbar_format}\""                                         >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "########### END USER INPUTS ###########"                                      >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo ""                                                                             >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
echo "set pm3d map"                                                                 >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
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
    line=`cat ${Input_Data_Folder}/output_filename_${Variables[$i]}_minmax.txt`
    linearr=( ${line} )
    minval=${linearr[0]}
    maxval=${linearr[1]}

    maxval=`expr $maxval*$colorbar_range_expander*${Variables_Normalizers[$i]}`
    minval=`expr $minval*$colorbar_range_expander*${Variables_Normalizers[$i]}`

    # Next line checks if ${Variables_Normalizers[$i]} is negative. If true, then Normalizer_sign=1.
    # Then we reverse min and max values on y-scale.
    Normalizer_sign=`echo "${Variables_Normalizers[$i]} < 0" | bc -l`
    if [ $Normalizer_sign -eq 1 ];then
       temp_minval=${maxval}
       temp_maxval=${minval}
       minval=${temp_minval}
       maxval=${temp_maxval}
    fi

    # label/enumerate the plots
    plot_number=$(( $i + 1 ))

    # generate the separate plots in gnuplot script
    if [ $(( $i + 1 )) -eq ${Nplots} ]; then
       echo "unset label"                                                           >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set label \"($plot_number)\" at -1.8,${SinglePlotHeight}/(${Nplus}+1)" >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "# Last plot with x-axis notation - ### Adjust notation here ###"       >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo ""                                                                      >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "# Set the x-axis notation for the bottom/last plot"                    >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "unset format"                                                          >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "unset xlabel"                                                          >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set xlabel \"${xlabel}\""                                              >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo ""                                                                      >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set ylabel \"${ylabel}\""                                              >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set cblabel \"${Variables[$i]}\""                                      >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set cbrange [${minval}:${maxval}]"                                     >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set xrange [${xstart}:${xstop}]"                                       >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set size SinglePlotSizeW,SinglePlotSizeH"                              >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set xtics ${xtics}"                                                    >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "splot \"${Input_Data_Folder}/output_filename_${Variables[$i]}.txt\" matrix using (\$1):(1+\$2):(${Variables_Normalizers[$i]}*\$3) notitle" >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo ""                                                                      >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
    else
       echo "unset label"                                                           >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set label \"($plot_number)\" at -1.8,${SinglePlotHeight}/(${Nplus}+1)" >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set ylabel \"${ylabel}\""                                              >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set cblabel \"${Variables[$i]}\""                                      >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set cbrange [${minval}:${maxval}]"                                     >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set xrange [${xstart}:${xstop}]"                                       >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set size SinglePlotSizeW,SinglePlotSizeH"                              >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "set xtics ${xtics}"                                                    >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
       echo "splot \"${Input_Data_Folder}/output_filename_${Variables[$i]}.txt\" matrix using (\$1):(1+\$2):(${Variables_Normalizers[$i]}*\$3) notitle" >> ${Output_Plot_Folder}/${Gnuplot_Script_name}
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
