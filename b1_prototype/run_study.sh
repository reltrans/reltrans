#!/bin/bash
export HEADAS=$(python3 -c "import xspectrampoline_helpers as h; print(h.get_HEADAS())")
cd "$(dirname "$0")"
run() { for g in "$@"; do for f in dc lo hi; do python3 study.py $g $f 2>&1 | grep -E "^$g |Error|Traceback" ; done; done; }
run default a05h3 a09i50h20 > study_A.log 2>&1 &
run i70 i10h10r3 > study_B.log 2>&1 &
wait
