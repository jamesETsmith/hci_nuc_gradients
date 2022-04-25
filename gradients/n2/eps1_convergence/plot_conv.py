import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append("../../../geometry_opt/hciscf/analysis/")
from plotting_utils import set_context, set_palette

os.makedirs("_figs", exist_ok=True)


def rms_error(g1: np.ndarray, g2: np.ndarray) -> float:
    return np.std((g1 - g2).flatten())


#
# Load Data
#
epsilon1 = [5e-3, 1e-3, 1e-4, 1e-5, 1e-6]
casscf = np.load("../single_point/_data/grad_casscf.npy")

vhciscf = [
    np.loadtxt("_data/vhciscf_{}.txt".format(eps1))
    for eps1 in epsilon1
]
vhciscf_err = [rms_error(shci, casscf) for shci in vhciscf]

vhciscf_aa = [np.loadtxt("_data/vhciscf_aa_{}.txt".format(eps1)) for eps1 in epsilon1]
vhciscf_aa_err = [rms_error(shci, casscf) for shci in vhciscf_aa]

hciscf = [
    np.loadtxt("_data/hciscf_{}.txt".format(eps1))
    for eps1 in epsilon1
]
hciscf_err = [rms_error(shci, casscf) for shci in hciscf]

hciscf_aa = [np.loadtxt("_data/hciscf_aa_{}.txt".format(eps1)) for eps1 in epsilon1]
hciscf_aa_err = [rms_error(shci, casscf) for shci in hciscf_aa]


# Debug print
print("CASSCF Gradient Norm", np.linalg.norm(casscf))
print(vhciscf_err)
print(vhciscf_aa_err)
print(hciscf_err)
print(hciscf_aa_err)

#
# Plot
#
label = ["vHCISCF", "vHCISCF + AA", "HCISCF", "HCISCF + AA"]

plt.figure()
# sns.set_palette(["#" + x for x in ["ff595e", "ffca3a", "1982c4", "6a4c93"]])
# sns.set_style("ticks")
set_context("paper", font_scale=1.25)
set_palette(4)

# Plot Gaussian tolerances
gau_tol = [1.7e-3, 3e-4, 1e-5]
gau_label = ["Gau. Loose", "Gau. Default", "Gau. Tight"]
styles = ["-", "--", "-."]

for i, tol in enumerate(gau_tol):
    plt.gca().axhline(
        tol, label=gau_label[i], linewidth=2, linestyle=styles[i], color="k"
    )

plt.loglog(epsilon1, vhciscf_err, "o-", label=label[0])
plt.loglog(epsilon1, vhciscf_aa_err, "o-", label=label[1])
plt.loglog(epsilon1, hciscf_err, "o-", label=label[2])
plt.loglog(epsilon1, hciscf_aa_err, "o-", label=label[3])


# plt.xlim((epsilon1[0] * 1.1, epsilon1[-1]))
plt.xlim((6e-3, 9e-7))
plt.xlabel(r"$\epsilon_1$ (Ha)")
plt.ylabel("RMS Gradient Error (Ha/Bohr)")
plt.title(r"Gradient Error as a Function of $\epsilon_1$")
# plt.ylim((1e-6, 3e-2))
plt.legend(framealpha=1.0)

plt.tight_layout()
plt.savefig("_figs/n2_eps1_conv.pdf")
plt.savefig("_figs/n2_eps1_conv.png", dpi=600)
