# List MOD14A1/MYD14A1 C6.1 granules (LP DAAC cloud) for the 26 MODIS tiles of the cloud-free analysis, 2002-12-01..2022-12-31, from NASA CMR (no login needed).
# Usage: python cmr_list.py tiles.csv urls.csv
import sys, json, urllib.request, urllib.parse, csv, time
tiles = list(csv.DictReader(open(sys.argv[1]))); out = open(sys.argv[2], 'w'); out.write('sensor,h,v,name,url\n')
n = 0
for t in tiles:
    h, v = int(t['h']), int(t['v'])
    for prod in ['MOD14A1', 'MYD14A1']:
        q = {'short_name': prod, 'version': '061', 'temporal': '2002-12-01T00:00:00Z,2022-12-31T23:59:59Z', 'page_size': 2000,
             'readable_granule_name': f'{prod}.A*.h{h:02d}v{v:02d}.061*', 'options[readable_granule_name][pattern]': 'true'}
        url = 'https://cmr.earthdata.nasa.gov/search/granules.json?' + urllib.parse.urlencode(q)
        for k in range(5):
            try:
                d = json.load(urllib.request.urlopen(url, timeout=120)); break
            except Exception as e:
                print('retry', prod, h, v, e, flush=True); time.sleep(10 * (k + 1))
        names = set()
        for e in d['feed']['entry']:
            hrefs = [l['href'] for l in e.get('links', []) if l['href'].startswith('https') and l['href'].endswith('.hdf')]
            if not hrefs or e['title'] in names: continue
            names.add(e['title']); out.write(f"{prod},{h},{v},{e['title']},{hrefs[0]}\n"); n += 1
        print(prod, h, v, len(names), flush=True)
out.close(); print('total', n)
