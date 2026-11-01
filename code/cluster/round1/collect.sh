#!/bin/bash
# pack the results for copying to the local working folder
B=${RFNC_BASE:-/scratch/$USER/rfnc}; cd $B
tar czf $B/rfnc_results_cmip6.tgz -C $B/out cmip6
cp $B/work/fire/fire_*.npz $B/work/rain/rain_*.npz $B/ 2>/dev/null
tail -n 3 logs/*.log > $B/logs_tail.txt
ls -la $B/*.tgz $B/*.npz
