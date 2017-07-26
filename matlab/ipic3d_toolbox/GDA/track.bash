#!/bin/bash
echo "inizio track" |mail -s "from mmsids"  valsusa@gmail.com

unset DISPLAY

echo "starting matlab"
#matlab 
matlab > matlab_movie.out  2>&1 << EOF

addpath '/home/gianni/matlab3'

track
exit

EOF

echo "fine track" |mail -s "from mmsids"  valsusa@gmail.com
