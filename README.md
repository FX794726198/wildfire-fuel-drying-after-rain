# Wildfire fuels dry out faster after rain as the air over land warms and dries

Code and minimal derived data for the manuscript "Wildfire fuels dry out faster after rain as the air over land warms and dries" by Xiao Feng, Luigi Lombardo, Juan Du and Bo Chai (submitted to *Nature Communications*).

The repository holds every analysis and figure script and the derived data needed to reproduce the figures and tables. The input data are public (see Data availability in the manuscript) and are not redistributed here. The scripts that download or extract them are included. Large intermediate files are left out because they can be rebuilt from the public data with the scripts (see [Not included](#not-included)).

## Quick start

```bash
git clone https://github.com/FX794726198/wildfire-fuel-drying-after-rain
cd wildfire-fuel-drying-after-rain
python figs/fig1.py          # also fig2.py, fig3.py, fig4.py
python figs/figS.py          # Figs S1-S4 and S6; figS5st.py draws Fig. S5, figS67.py Figs S7 and S8
python figs/tables.py        # Tables S1 and S2 (base rows); tableS4.py and tableS5.py write Tables S4 and S5
```

Run every script from the repository root. The figures are written to `figs/` and the tables to `tables_out/` (set `PARTS` to change this). All twelve figures redrawn from `figs/data/` in a fresh copy of this repository are identical, pixel for pixel, to those in the manuscript.

Python 3 with numpy, pandas, scipy, statsmodels, matplotlib and cartopy is enough for the figures, tables and the cloud-free fire analysis. The full chain also needs xarray, zarr, gcsfs, cftime, lightgbm and h5py; `cloudfree/wc_cover.py` needs rasterio, and the fire-mask extraction needs pyhdf and curl.

## Layout

| Folder or file | What it holds |
|---|---|
| `code/` | Analysis scripts, from data extraction to the summary tables, and `cells_sel.csv` (the 1,704 analysis cells), `smon_obs_sel.npz` (observed fire-season months per cell) and `fs_jobs.csv` (the 21 CMIP6 runs). |
| `code/cluster/` | Batch scripts for the DelftBlue cluster: CMIP6 runs, fire detections and rain products (`round1/`), the CMIP6 summary step (`round2/`) and the Nelson dead fuel moisture model (`nelson/`). |
| `cells_all.csv`, `cells_sel.csv` | All 1.5° land cells and the analysis cells. |
| `out_era5/` | Per-cell results for the main analysis and every alternative input (`cells_<variant>.csv`), and the drying time after each of the 144,602 rain events of the main analysis (`events_obs.csv.gz`) and of the CHIRPS version (`events_rain_chirps.csv.gz`). |
| `res/` | Summary tables of all results quoted in the manuscript, including the station results (`isd_*`), the trend-removal and DAMIP results (`fixed_*`), the climate-model comparison (`cmip6_*`, `samecov*`, `models_*`), the fire-detection analysis (`fire_*`) and the Nelson model (`nelson_*`). |
| `cloudfree/` | Fire activity relative to cloud-free observations after rain in five regions (Fig. 2e–g, Fig. S3c–e,g–l, Section S3): scripts, the daily fire-mask counts per cell (`fm.npz`), tree cover per cell (`cell_cover.csv`), ERA5 24-hour rain and daily FFMC for the 347 cells (`era5g_p24/`, `ffmc_p347.npz`), the MODIS fire-mask extraction (`extraction/`) and the results (`res/`). |
| `figs/` | Figure and table scripts; `figs/data/` holds the data plotted in every figure. |
| `tables/` | Supplementary Tables S1–S5 as in the manuscript (Markdown). |

## Where each result comes from

| Manuscript item | Scripts (in run order) | Output |
|---|---|---|
| Drying time after rain from ERA5 (Fig. 1, Table S1) | `e1_era5_extract.py`, `e2_ffmc_global.py` (functions), `e2_run.py`, `a1_core.py`, then `figs/p1_data.py`, `figs/fig1.py`, `figs/tables.py` | `out_era5/cells_obs.csv`, `out_era5/events_obs.csv.gz`, `res/series_global.csv`, `figs/data/fig1_*` |
| Other rain products, wind and thresholds (Fig. 1b, Fig. S2, Table S1) | `e5_rain.py`, `e2_run.py` with the alternative inputs | `out_era5/cells_rain_chirps.csv`, `cells_rain_cpc.csv`, `cells_windobs.csv`, `cells_rt90.csv`, `cells_txrn.csv` |
| Nelson fuel-stick model (Fig. 1b, Fig. S2e,f) | `code/cluster/nelson/`, then `a6_nelson.py` | `res/nelson_events.csv.gz`, `res/nelson_summary.csv` |
| Fire detections (Fig. 2a–d, Fig. S3a,b, Section S3) | `e4_fire.py`, `a4_fire.py`, `a4b_ci.py`, `a4c_regci.py`, then `figs/p2_data.py`, `figs/fig2.py` | `res/fire_*`, `figs/data/fig2_*` |
| Gradient-boosting check (Fig. S3f) | `a4_fire.py`, `a4b_ci.py` | `res/fire_summary_modis.json`, `figs/data/fig2_ci_modis.json` |
| Fire activity relative to cloud-free observations (Fig. 2e–g, Fig. S3c–e,g–l, Section S3) | `cloudfree/extraction/`, `cloudfree/combine_fm.py`, `cloudfree/check_day_alignment.py`, `cloudfree/wc_cover.py`, then `cloudfree/cloudfree_analysis.py <group>` and `cloudfree/cloudfree_analysis.py merge`; then `figs/fig2.py` and `figs/figS.py S3` | `cloudfree/fm.npz`, `cloudfree/cell_cover.csv`, `cloudfree/res/*`, `figs/data/cloudfree_results.json`, `figs/data/figS3_cloudfree_cells.csv` |
| Removing temperature and humidity trends after the same events (Fig. 3b, Fig. S4, Table S2) | `e7_fixed.py`, `a11_fixed_summary.py`, then `figs/p3_data.py`, `figs/ps_data.py`, `figs/p3b_fixed.py` (in this order) | `res/fixed_noon_cap31.csv`, `res/fixed_chirps_cap31.csv`, `figs/data/fig3_decomp.csv`, `figs/data/figS3_maps.csv` |
| Drying relation (Fig. 3c) | `a2_scaling.py`, `figs/p3_data.py`, `figs/fig3.py` | `res/scaling_cells.csv`, `res/scaling_regions.csv`, `figs/data/fig3_*` |
| Fire-season days with dry fuels (Methods, Discussion) | `a7_decomp.py`, `a7b_dryday_attrib.py` | `res/decomp_obs.*`, `res/dryday_attrib*` |
| Station checks with noon values (Fig. S5a–c) | `a8_isd_list.py`, `a8_isd_select2.py`, `a8_isd_get.py`, `e1s_era5_station_hours.py`, `a8_isd_ffmc.py`, `a8b_isd_summary.py`, then `figs/ps5_data.py`, `figs/figS5st.py` | `res/isd_main*`, `res/isd_cpc*`, `figs/data/figS5_*` |
| Station checks with daily maximum temperature and minimum humidity, versions A and B (Fig. 1b, Fig. S5d, Table S1, Section S5) | `a10_isd_get3.py`, `e1t_era5_txrn3.py`, `a12_isd_txrn.py`, then `figs/ps5b_data.py` | `res/isd_txrn.csv`, `figs/data/figS5_txrn.json` |
| Human influence from DAMIP (Fig. 4a, Fig. S6) | `e6b_damip.py` (`run_e6b.sh`), `e1m_era5_dailymean.py`, `e7_fixed.py`, `a11_fixed_summary.py`, `a14_damip_rh.py`, then `figs/p4c_data.py`, `figs/figS.py` | `res/fixed_noon_cap31_dmp.csv`, `res/fixed_noon_cap31_to2020_dmp.csv`, `res/damip_rh.csv`, `figs/data/fig4_damip.csv`, `figs/data/figS6_damip_TH.csv` |
| CMIP6 with the observed fire season in the same cells (Fig. 4d,e, Table S5, Section S8) | `e3_cmip6.py` (`run_fs.sh`, `fs_jobs.csv`), then `final_models.sh` and `a13c_diag.py` | `res/samecov*`, `res/models_table.*`, `res/models_diag.json`, `figs/data/fig4_samecov*`, `figs/data/fig4_obs_ci.json` |
| CMIP6 with each model's own fire season, control runs and warming levels (Fig. 4b,c,e, Figs S7 and S8, Table S4) | `e3_cmip6.py` (`code/cluster/round1/job_cmip6.sh`), `a5_cmip6.py`, then `figs/p4_data.py`, `figs/tableS4.py`, `figs/figS67.py` | `res/cmip6_*`, `figs/data/fig4_*`, `figs/data/figS6_series.csv` |
| Large ensembles (Fig. 4b) | `run_le2.sh`, `a9_le.py` | `res/le_members.csv`, `figs/data/fig4_le_members.csv` |
| Supplementary figures S1–S4 | `figs/ps_data.py`, `figs/figS.py` | `figs/data/figS1_*` to `figS4_*` |

`figs/tables.py` writes the base rows of Tables S1 and S2. The station rows of Table S1 with daily maximum temperature and minimum humidity are taken from `figs/data/figS5_txrn.json`, and the paired differences of Table S2 from `res/fixed_noon_cap31.csv` and `res/fixed_chirps_cap31.csv`. The two rows of the GHCN station check (510 stations with station rain) are read from `res/robustness_versions.csv`. The final tables are in `tables/`.

The global warming levels of the CMIP6 models are taken from the public table of M. Hauser (https://github.com/mathause/cmip_warming_levels, file `cmip6_warming_levels_all_ens_1850_1900.csv`, reference period 1850–1900). `a5_cmip6.py` reads it from `gwl/` in the folder where it is run.

## Not included

These intermediate files are large and are rebuilt from public data by the scripts named here.

| Folder or file | Contents | Script |
|---|---|---|
| `era5g/` | Daily ERA5 noon weather, 24-hour rain, daily maximum temperature and minimum humidity | `e1_era5_extract.py` |
| `era5m/` | Daily mean ERA5 temperature and humidity | `e1m_era5_dailymean.py` |
| `out_era5/ffmc_obs.npz`, other `out_era5/events_*.csv.gz` | Daily FFMC from ERA5 and the event lists of the alternative inputs | `e2_run.py` |
| `fx/` | Drying times after the same rain events with modified weather | `e7_fixed.py` |
| `fs/` | Output of the 21 CMIP6 models run with the observed fire-season months | `e3_cmip6.py` (`run_fs.sh`) |
| `dbo/` | Gridded MODIS and VIIRS fire detections, GFED emissions, CHIRPS and CPC rain | `e4_fire.py`, `e5_rain.py` |
| `isd/` | ISD station series and ERA5 at the stations | `a8_isd_*.py`, `a10_isd_get3.py`, `e1s_era5_station_hours.py`, `e1t_era5_txrn3.py` |
| `damip2/`, `le/` | DAMIP forced changes and large-ensemble members | `e6b_damip.py` (`run_e6b.sh`), `run_le2.sh` |
| MODIS fire masks | Daily MOD14A1 and MYD14A1 granules, about 46,000 files | `cloudfree/extraction/` |

## The cloud-free fire-activity analysis (`cloudfree/`)

- Cells: `cells_cloudfree.csv` lists the 347 cells of the five regions (SAM, South American monsoon; SES, southeastern South America; ESAF, eastern southern Africa; NAU, northern Australia; CAF, central Africa).
- Fire masks: `extraction/cmr_list.py tiles.csv urls.csv` lists the MOD14A1 and MYD14A1 Collection 6.1 granules of the 26 MODIS tiles from NASA CMR. `extraction/extract_fm.py` downloads each granule with curl, using your own NASA Earthdata login in `~/.netrc`, and counts the pixels of each cell, day and satellite by fire-mask class; `job_fm.slurm` is an example batch job. `combine_fm.py` combines the outputs into `fm.npz`: fire pixels of nominal or high confidence (`f89`) and of all confidence levels (`f79`), cloud-free land pixels (`valid`) and cloudy land pixels (`cloud`), each as [cell, day from 1 December 2002, satellite].
- ERA5 rain and FFMC: the analysis reads the 24-hour rain from `era5g/era5_<year>.npz` and the FFMC from `out_era5/ffmc_obs.npz`. To run it without rebuilding these, link the included rain subset first (`ln -s cloudfree/era5g_p24 era5g`); when `out_era5/ffmc_obs.npz` is absent, the script reads the FFMC of the 347 cells from `cloudfree/ffmc_p347.npz`.
- Analysis: `python cloudfree/cloudfree_analysis.py <group> [--boot 1000]` writes `cloudfree/res/<group>.json`. The groups are `pooled`, `regions`, `periods`, `leave_one_out`, `conditions` and `cover`; `all` runs all 24 analyses. `cloudfree_analysis.py merge` collects them into `res/cloudfree_results.json` and `res/cloudfree_summary.csv`. With 1,000 resamples the pooled analyses take about 30–90 minutes each on one core.
- `fastglm.py` fits the Poisson regression by iteratively reweighted least squares in numpy, with the same estimates as statsmodels and much faster in the resampling; `FASTGLM=0` uses statsmodels instead.

## Definitions

- Drying time: days from a fire-season rain of at least 10 mm, falling when the FFMC of the previous day was at least 85, until the FFMC first reaches 85 again (1–31 days). The model comparison caps drying times at 15 days and pools the events of each cell and period.
- Intervals come from resampling 5° blocks of cells or stations (`figs/figcommon.py`, `block_ci`), and in the cloud-free analysis from resampling blocks of 5° and one year.
- DAMIP responses: anthropogenic = historical (with SSP2-4.5 after 2014) minus hist-nat; greenhouse gases and aerosols = the single-forcing run minus its own 1850–1900 mean.
- Station version A screens the station values against ERA5 and fills gaps with ERA5. Version B uses only the station's own records.
- Changes of relative humidity are given in % (units of relative humidity).

## Licence

Code: MIT licence (see `LICENSE`). Derived data: Creative Commons Attribution 4.0 International (CC BY 4.0). Please cite the paper when using either.
