import sys
import numpy as np
import rmsd

sys.path.append("..")
from analysis_utils import parse_mcscf_energies, parse_multi_XYZ

# Check energies as a sanity check
# fmt: off
natorb = parse_mcscf_energies("../../dz/1_A/40e_40o_eps1_response/._logs_bkp/_logs/opt_5e-05.out")
no_natorb = parse_mcscf_energies("../../dz/1_A/40e_40o_eps1_response/_logs/opt_5e-05.out")
max_diff = 1e-14
for i in range(min(natorb.size,no_natorb.size)):
    diff = abs(no_natorb[i]-natorb[i])
    max_diff = max(max_diff, diff)
    print(f"w/o NO: {no_natorb[i]:.6f}\tw/ NO: {natorb[i]:.6f}\t Diff: {diff:.1e}")
# fmt: on
if max_diff > 1e-6:
    print(f"WARNING: BIG DIFFERENCE IN ENERGIES ({max_diff:.1e})")

# Gradients

# fmt: off
# TODO THIS IS THE BEST ANSWER WITH NATORBS
best_answer = np.loadtxt("../../dz/1_A/40e_40o_eps1_response/._data_bkp/_data/1_A_DZ_40e_40o_eps1=3e-05_gradients.txt")
no_natorb = np.loadtxt("../../dz/1_A/40e_40o_eps1_response/_data/1_A_DZ_40e_40o_eps1=5e-05_gradients.txt")
natorb = np.loadtxt("../../dz/1_A/40e_40o_eps1_response/._data_bkp/_data/1_A_DZ_40e_40o_eps1=5e-05_gradients.txt")
# fmt: on

print(f"Effect of Natural Orbitals on {np.linalg.norm(no_natorb - natorb):.1e}")
print(f"Difference from best answer w/  {np.linalg.norm(best_answer - natorb):.1e}")
print(f"Difference from best answer w/o {np.linalg.norm(best_answer - no_natorb):.1e}")

# Compare Gradients without Natorb
# fmt: off
eps1_3e_5 = np.loadtxt("../../dz/1_A/40e_40o_eps1_response/_data/1_A_DZ_40e_40o_eps1=3e-05_gradients.txt")
eps1_5e_5 = np.loadtxt("../../dz/1_A/40e_40o_eps1_response/_data/1_A_DZ_40e_40o_eps1=5e-05_gradients.txt")
eps1_7_5e_5 = np.loadtxt("../../dz/1_A/40e_40o_eps1_response/_data/1_A_DZ_40e_40o_eps1=7.5e-05_gradients.txt")
eps1_1e_4 = np.loadtxt("../../dz/1_A/40e_40o_eps1_response/_data/1_A_DZ_40e_40o_eps1=0.0001_gradients.txt")
eps1_2e_4 = np.loadtxt("../../dz/1_A/40e_40o_eps1_response/_data/1_A_DZ_40e_40o_eps1=0.0002_gradients.txt")
# fmt: on

n_atom = eps1_5e_5.shape[0]
# fmt: off
print(f"Difference in 3e-5 and 5e-5   gradients {np.linalg.norm(eps1_3e_5-eps1_5e_5)/ np.sqrt(n_atom)}")
print(f"Difference in 3e-5 and 7.5e-5 gradients {np.linalg.norm(eps1_3e_5-eps1_7_5e_5)/ np.sqrt(n_atom)}")
print(f"Difference in 3e-5 and 1e-4   gradients {np.linalg.norm(eps1_3e_5-eps1_1e_4)/ np.sqrt(n_atom)}")
print(f"Difference in 3e-5 and 2e-4   gradients {np.linalg.norm(eps1_3e_5-eps1_2e_4)/ np.sqrt(n_atom)}")
# fmt: on

exit(0)

# Compare geometries from tight geom opt
state = "3_B"
_, geoms_tight, _ = parse_multi_XYZ(
    f"../../dz/{state}/40e_40o_tight/_data/{state}_DZ_40e_40o_eps1=5e-05.xyz"
)

_, geoms_loose, _ = parse_multi_XYZ(
    f"../../dz/{state}/40e_40o/_data/{state}_DZ_40e_40o_eps1=7.5e-05.xyz"
)

# rmsd_factor = np.sqrt(geoms_tight[0].shape[0])

for i in range(min(len(geoms_loose), len(geoms_tight))):
    # diff = np.linalg.norm(geoms_tight[i] - geoms_loose[i]) / rmsd_factor
    diff = rmsd.kabsch_rmsd(geoms_tight[i], geoms_loose[i])
    print(f"Difference in geom at iter {i:2d} = {diff:.2e}")

#
diff = rmsd.kabsch_rmsd(geoms_tight[-1], geoms_loose[-1])
print(f"Difference in geom at iter FINAL = {diff:.2e}")
