#!/bin/sh
Xvfb :7

echo "starting matlab"

export DISPLAY=:7

matlab  > Etot3Dfluxsweep.out <<EOF
energia
ohm
EOF
