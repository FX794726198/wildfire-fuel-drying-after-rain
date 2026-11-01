#!/bin/bash
# Build the Nelson dead fuel moisture model from NFDRS4 (public domain, https://github.com/firelab/NFDRS4) as a small shared library.
set -e
B=${RFNC_BASE:-/scratch/$USER/rfnc}; cd $B
[ -d NFDRS4 ] || git clone --depth 1 https://github.com/firelab/NFDRS4.git
S=$B/NFDRS4/lib/NFDRS4
cp $B/code2/dfm_c.cpp $S/src/
g++ -O2 -shared -fPIC -I$S/include -o $B/code2/libdfm.so $S/src/dfm_c.cpp $S/src/deadfuelmoisture.cpp $S/src/dfmcalcstate.cpp
ls -la $B/code2/libdfm.so
