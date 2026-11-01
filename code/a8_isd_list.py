# List isd-lite station-year files (key, size) for 1979-2022 from the public NOAA S3 bucket.
import requests, re, sys, os
from concurrent.futures import ThreadPoolExecutor
B = 'https://noaa-isd-pds.s3.amazonaws.com/'; out = './isd/lists'
def year(y):
    f = f'{out}/{y}.tsv'
    if os.path.exists(f): return
    rows, tok = [], None
    while True:
        p = {'list-type': '2', 'prefix': f'isd-lite/data/{y}/', 'max-keys': '1000'}
        if tok: p['continuation-token'] = tok
        for k in range(6):
            try: r = requests.get(B, params=p, timeout=60); r.raise_for_status(); break
            except Exception: pass
        x = r.text
        rows += re.findall(r'<Key>isd-lite/data/\d+/([^<]+)\.gz</Key>.*?<Size>(\d+)</Size>', x)
        m = re.search(r'<NextContinuationToken>([^<]+)</NextContinuationToken>', x)
        if not m: break
        tok = m.group(1)
    open(f, 'w').write(''.join(f'{a}\t{b}\n' for a, b in rows)); print(y, len(rows), flush=True)
with ThreadPoolExecutor(11) as ex: list(ex.map(year, range(1979, 2023)))
