#!/bin/bash
echo "inizio sat traces" |mail -s "from mmsids on $OF"  valsusa@gmail.com

unset DISPLAY

echo "starting matlab"
#matlab 
matlab > matlab_movie.out  2>&1 << EOF

addpath '/home/gianni/matlab3'

movie_ve
exit

EOF

convert -density 300 Ven_combo*.eps Ven_combo.gif
convert -density 300 Vel_combo*.eps Vel_combo.gif
convert -density 300 Vem_combo*.eps Vem_combo.gif
convert -density 300 Ven0*.eps Ven.gif
convert -density 300 Vel0*.eps Vel.gif
convert -density 300 Vem0*.eps Vem.gif
