#!/bin/bash
echo "inizio sat traces" |mail -s "from mmsids on $OF"  valsusa@gmail.com
unset DISPLAY
echo "starting matlab"

matlab  > matto.out 2>&1 << EOF
tred54
tred60
EOF


echo "fine sat traces" |mail -s "from mmsids on $OF"  valsusa@gmail.com
