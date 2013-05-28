reset
# Sum and plot together two indentical matrices (must be equal row, column)
############# USER INPUTS #############

# Number of stacked plots
Nplots=1
SinglePlotHeight=0.8 #800
plotwidth=1 #1000

# Define the total plot size
plotheight=Nplots*SinglePlotHeight
set term postscript eps enhanced #set terminal jpeg size plotwidth,plotheight
set output "./test.eps"

# In order to leave room for axis and tic labels underneath, a N+1 or N+2 -plot layout is created but only use the top N slots.

Nplus=Nplots+1
set multiplot layout Nplus,1 title "log(|FFT(Ey)|) along X (96 probes->48 FFT points)" font "Helvetica,20"
set tmargin 0
set bmargin 0

leftmargin=5
rightmargin=5

time_axis_scaler=0.0097*0.125 #Bx0*dt # the 2D data is written according to cycle number and this is to convert it to [wci*t] units
space_axis_scaler=1

coord11=00.0
coord12=16771.0
coord21=00.0
coord22=48.0

delta_X_between_probes=0.3906  #in units d_i : tred46 we have dx=0.3906 ,dy=0.46875, dz=0.39063
ks=2.0*3.1415927/(48.0*delta_X_between_probes)

scaleX=time_axis_scaler
scaleY=ks

########### END USER INPUTS ###########

set pm3d map

# All plots except the bottom one will not have x-axis notation
unset format

## Start form the top - First Figure

set xlabel "w_{ci}t"
set ylabel "k_m/d_i"
set cbrange [0.000001:0.0015]
set palette model RGB # defined (0 "white", 1 "black")
set palette defined
set logscale zcb
set format cb "%L"
set cblabel "log(|FFT(Ey(x))|)"
set xrange [coord11*scaleX:20] #coord12*scaleX]
set yrange [coord21*scaleY:coord22*scaleY]
splot "fft_B_tred46_along_X.txt" matrix using (scaleX*$1):(scaleY*$2):(1.0*$3) notitle

unset multiplot
reset
