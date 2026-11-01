# Check one downloaded MOD14A1/MYD14A1 file: SDS names, shapes, and FireMask class x QA land/water state (bits 0-1) cross-table.
# Usage: python qa_check.py file.hdf
import sys, numpy as np
from pyhdf.SD import SD, SDC
s = SD(sys.argv[1], SDC.READ); print('SDS:', {k: v[1] for k, v in s.datasets().items()})
fm = s.select('FireMask')[:]; qa = s.select('QA')[:]; print('FireMask', fm.shape, fm.dtype, 'QA', qa.shape, qa.dtype)
for k in range(min(2, fm.shape[0])):
    t = np.zeros((10, 4), int)
    np.add.at(t, (fm[k].ravel().clip(0, 9), (qa[k].ravel() & 3)), 1)
    print(f'layer {k}: rows FireMask 0-9, cols QA&3 = 0 water,1 coast,2 land,3 unused'); print(t)
    print('QA bit2 (day flag) counts:', np.bincount(((qa[k] >> 2) & 1).ravel(), minlength=2))
