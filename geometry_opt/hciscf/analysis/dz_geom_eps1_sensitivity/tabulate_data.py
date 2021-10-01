import sys
import os
import numpy as np
import pandas as pd
from rmsd import kabsch_rmsd

sys.path.append("..")
from analysis_utils import parse_multi_XYZ

os.makedirs("_data", exist_ok=True)

# fmt: off
eps_loose = parse_multi_XYZ("../../dz/1_A/40e_40o/_data/1_A_DZ_40e_40o_eps1=7.5e-05.xyz")
eps_tight = parse_multi_XYZ("../../dz/1_A/40e_40o_tight/_data/1_A_DZ_40e_40o_eps1=5e-05.xyz")
# fmt: on

min_iter = min(len(eps_loose[0]), len(eps_tight[0]))
for i in range(min_iter):
    diff = kabsch_rmsd(eps_loose[1][i], eps_tight[1][i])
    print(f"Iter. {i:2d} RMSD {diff:.1e} \u212B")


print("Final Geom RMSD")
print(kabsch_rmsd(eps_loose[1][-1], eps_tight[1][-1]))
