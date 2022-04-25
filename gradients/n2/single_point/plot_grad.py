import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append("../../../geometry_opt/hciscf/analysis/")
from plotting_utils import set_context, set_palette


def rms_error(g1: np.ndarray, g2: np.ndarray) -> float:
    return np.std((g1-g2).flatten())


#
# Read Data
#
eps1 = 5e-3
casscf = np.load("_data/grad_casscf.npy")


vhciscf = np.load(f"_data/grad_vhciscf.npy")
vhciscf_aa = np.load(f"_data/grad_vhciscf_aa.npy")
hciscf = np.load(f"_data/grad_hciscf.npy")
hciscf_aa = np.load(f"_data/grad_hciscf_aa.npy")

n_atom = casscf.shape[0]
shci_grads = [vhciscf, vhciscf_aa, hciscf, hciscf_aa]
shci_errors = [rms_error(gr,casscf) for gr in shci_grads]

print("CASSCF Grad Norm", np.linalg.norm(casscf))
print("HCISCF Gradient Errors\n", shci_errors)
# for gr in shci_grads:
#     print(gr[0][2])


#
# Plot
#
tick_label = ["vHCISCF", "vHCISCF + AA", "HCISCF", "HCISCF + AA"]
os.makedirs("_figs", exist_ok=True)

plt.figure()

set_palette(4)
set_context("paper", font_scale=1.25)

sns.barplot(np.arange(len(shci_errors)), shci_errors)

plt.xticks(np.arange(len(tick_label)), tick_label)

plt.ylabel("RMS Gradient Error (Ha/Bohr)")
plt.title("Gradient Error as a Function of Method Variant")
plt.tight_layout()
plt.savefig(f"_figs/n2_single_point_error.png", dpi=600)
plt.savefig(f"_figs/n2_single_point_error.pdf")
plt.close()
