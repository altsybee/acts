#!/bin/bash

python material_budget_gdml_scan.py \
  --gdml ../o2sim_geometry.gdml \
  --plot-vs eta --make-th2\
  --phi-min 0 --phi-max 6.283185307 --n-phi 3 \
  --eta-min -2.8 --eta-max 2.8 --n-eta 161 \
  --rmin 0 --rmax 86 --zmin -300 --zmax 300  \
  --origin-x 0 --origin-y 0 --origin-z 0 \
  --top-materials 10 --max-steps 20000 \
  --collect-physvol barrel_1 \
  --table-etas "0.01, " \
  --table-phis "0.01, " \
  --table-out-prefix myrays --pt 1.0 --bz 2.0 --charge 1 --arc-step-mm 1.0\
  --out mat_budget.png \
  --root-out mat_budget.root

  # --table-etas "  2.1,  " \
  # --table-phis "0.01,  " \
  # --table-etas "0, 0.1, 0.01, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8," \
  # --table-phis "0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01" \

  # --table-etas 0.0 \
  # --table-phis 0.0 \
  # --table-out-prefix ray_traversal \
  # --segments-threshold 0.0001 \
  # --merge-gap-mm 1.0 \
  # --merge-barrel-name barrel_1

  # --table-etas "0,1.5,2.2" --table-phis "0,0,0" --table-out-prefix myrays \
  
  