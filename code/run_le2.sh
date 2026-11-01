#!/bin/bash
# Large-ensemble CMIP6 runs (historical + SSP2-4.5 to 2022); a lock directory per member prevents duplicate runs.
cd .
tail -n +2 code/le_members4.csv | while IFS=, read m mem; do echo "$m $mem"; done | \
xargs -P ${1:-3} -n 2 sh -c '[ -f le/$0_hist245_$1.csv.gz ] && exit 0; mkdir le/lock/$0_$1 2>/dev/null || exit 0; YMAX=2022 MEMBER=$1 CAT=pangeo-cmip6.csv python3 code/e3_cmip6.py $0 hist245 code/cells_sel.csv le >> le/log_$0.txt 2>&1'
