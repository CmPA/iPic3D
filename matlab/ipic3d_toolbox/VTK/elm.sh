#!/bin/sh
Xvfb :7

echo "starting matlab"

export DISPLAY=:7

matlab  > elm.out <<EOF
elm_energy_eq
ohm
EOF
