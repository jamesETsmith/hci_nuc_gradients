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
# eps1 = 1e-3
casscf = np.load("_data/grad_casscf.npy")

for eps1 in [5e-4]:
    vhciscf = np.load(f"_data/grad_vhciscf_{eps1:.1e}.npy")
    vhciscf_aa = np.load(f"_data/grad_vhciscf_aa_{eps1:.1e}.npy")
    hciscf = np.load(f"_data/grad_hciscf_{eps1:.1e}.npy")
    hciscf_aa = np.load(f"_data/grad_hciscf_aa_{eps1:.1e}.npy")
    # vhciscf_aa = np.zeros_like(casscf)
    # hciscf = np.zeros_like(casscf)
    # hciscf_aa = np.zeros_like(casscf)

    shci_grads = [vhciscf, vhciscf_aa, hciscf, hciscf_aa]
    shci_errors = [rms_error(gr,casscf) for gr in shci_grads]

    print("CASSCF Grad Norm",np.linalg.norm(casscf))
    print("HCISCF Gradient Errors\n", shci_errors)
    # for gr in shci_grads:
    #     print(gr[0][2])


    #
    # Plot
    #
    tick_label = ["vHCISCF", "vHCISCF + AA", "HCISCF", "HCISCF + AA"]
    os.makedirs("_figs", exist_ok=True)

    plt.figure()
    # sns.set_palette(["#" + x for x in ["ff595e", "ffca3a", "1982c4", "6a4c93"]])
    # sns.set_style("ticks")
    set_palette(4)
    set_context("paper", font_scale=1.25)

    sns.barplot(np.arange(len(shci_errors)), shci_errors)
    plt.xticks(np.arange(len(tick_label)), tick_label)

    plt.ylabel("RMS Gradient Error (Ha/Bohr)")
    plt.title("Gradient Error as a Function of Method Variant")
    plt.tight_layout()
    plt.savefig(f"_figs/ni_edt2_single_point_error_{eps1:.1e}.png", dpi=600)
    plt.savefig(f"_figs/ni_edt2_single_point_error_{eps1:.1e}.pdf")
    plt.close()
