#!/bin/bash
# v1.2: same-coverage comparison, per-model table, Fig. 4 data and figure once all 21 fixed-season runs exist
cd .
python3 code/a13_samecov.py && python3 code/a13b_models_table.py > res/models_table.txt && python3 figs/tableS5.py > /dev/null && python3 figs/p4c_data.py > /dev/null && cd figs && python3 fig4.py
