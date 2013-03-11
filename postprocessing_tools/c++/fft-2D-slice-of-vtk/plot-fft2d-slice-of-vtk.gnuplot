reset

############# USER INPUTS #############

# Number of stacked plots
Nplots=2
SinglePlotHeight=400
plotwidth=1000

# Define the total plot size
plotheight=Nplots*SinglePlotHeight
set term postscript eps enhanced
set output "./fft-2D-slice-fragment-VTK-xz.eps"

# In order to leave room for axis and tic labels underneath, a N+1 or N+2 -plot layout is created but only use the top N slots.

Nplus=Nplots+1
set multiplot layout Nplus,1 title "FFT of 2D slice fragment of rho_tot_electrons" font "Helvetica,20"
set tmargin 0
set bmargin 1

leftmargin=5
rightmargin=5

coord11=10.0
coord12=30.0
coord21=00.0
coord22=15.0

index11=135
index12=407
index21=0
index22=135

scaleX=(30.0-10.0)/(407-135)
scaleY=(15.0-00.0)/(135-0)

scaleXfft=2*3.1415927/(407-135)*0.125/0.0097
scaleYfft=2*3.1415927/((135-0)/2+1)*0.125/0.0097

########### END USER INPUTS ###########

set pm3d map

# All plots except the bottom one will not have x-axis notation
unset format

## Start form the top - First Figure

set xlabel "X "
set ylabel "Z "
set cbrange [:]
set palette color
set palette model RGB
set palette defined
set format cb "%g"
set xrange [135*scaleX:407*scaleX]
set yrange [0*scaleY:135*scaleY]
splot "./output_2D_slice_VTK.txt" matrix using (10.0+scaleX*$1):(00.0+scaleY*$2):(1.0*$3) notitle

set xlabel "w/wpi"
set ylabel "w/wpi"
set cbrange [:]
set palette color
set palette model RGB
set palette defined
set logscale zcb
set format cb "%L^{10}"
set xrange [0*scaleXfft:(index12-index11)*scaleXfft]
set yrange [0*scaleYfft:((index22-index21)/2+1)*scaleYfft]
splot "./output_fft_2D_slice_VTK.txt" matrix using (scaleXfft*$1):(scaleYfft*$2):(1.0*$3) notitle

unset multiplot
reset
