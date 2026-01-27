@!/bin/bash

make -f revmakefile clean

make -f revmakefile lib

read -p 'Press Enter '

# python3 test_single_model.py
python3 test_plunging_regions_plots.py
