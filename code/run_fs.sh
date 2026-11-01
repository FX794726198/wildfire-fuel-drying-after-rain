#!/bin/bash
# 21 CMIP6 models (first member, historical + SSP to 2022) with the observed fire-season months (code/smon_obs_sel.npz).
cd .
tail -n +2 code/fs_jobs.csv | tr ',' ' ' | \
xargs -P ${1:-4} -n 3 sh -c '[ -f fs/$0_$1_$2.csv.gz ] && exit 0; mkdir fs/lock/$0 2>/dev/null || exit 0; YMAX=2022 SMON=code/smon_obs_sel.npz MEMBER=$2 CAT=pangeo-cmip6.csv python3 code/e3_cmip6.py $0 $1 code/cells_sel.csv fs >> fs/log_$0.txt 2>&1'
