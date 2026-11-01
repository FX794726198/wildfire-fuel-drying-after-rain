cd .
for m in MRI-ESM2-0 CanESM5 CNRM-CM6-1 IPSL-CM6A-LR GFDL-ESM4 HadGEM3-GC31-LL CESM2; do echo $m; done | \
xargs -P 2 -I{} sh -c 'CAT=pangeo-cmip6.csv python3 code/e6b_damip.py {} cells_all.csv damip2 2>&1 | tail -2' > damip2/e6b.log 2>&1
echo DONE >> damip2/e6b.log
