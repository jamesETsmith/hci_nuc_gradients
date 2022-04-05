import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#
# Load Data
#
epsilon1 = [5e-3, 1e-3, 1e-4, 1e-5, 1e-6]
casscf = np.load("../single_point/_data/grad_casscf.npy")

vhciscf = [np.loadtxt("_data/vhciscf_{}.txt".format(eps1)) for eps1 in epsilon1]
vhciscf_err = [np.linalg.norm(casscf - shci) for shci in vhciscf]

# vhciscf_aa = [np.loadtxt("_data/vhciscf_aa_{}.txt".format(eps1)) for eps1 in epsilon1]
# vhciscf_aa_err = [np.linalg.norm(casscf - shci) for shci in vhciscf_aa]


hciscf = [np.loadtxt("_data/hciscf_{}.txt".format(eps1)) for eps1 in epsilon1]
hciscf_err = [np.linalg.norm(casscf - shci) for shci in hciscf]

# hciscf_aa = [np.loadtxt("_data/hciscf_aa_{}.txt".format(eps1)) for eps1 in epsilon1]
# hciscf_aa_err = [np.linalg.norm(casscf - shci) for shci in hciscf_aa]


# Debug print
# print(vhciscf_err)
# print(vhciscf_aa_err)
# print(hciscf_err)
# print(hciscf_aa_err)

#
# Plot
#
label = ["vHCISCF", "vHCISCF + AA", "HCISCF", "HCISCF + AA"]
# label = ["vHCISCF", "HCISCF",]

plt.figure()
sns.set_palette(["#" + x for x in ["ff595e", "ffca3a", "1982c4", "6a4c93"]])
sns.set_style("ticks")
sns.set_context("talk")

plt.loglog(epsilon1, vhciscf_err, "o-", label=label[0])
# plt.loglog(epsilon1, vhciscf_aa_err, "o-", label=label[1])
plt.loglog(epsilon1, hciscf_err, "o-", label=label[2])
# plt.loglog(epsilon1, hciscf_aa_err, "o-", label=label[3])


# Plot Gaussian tolerances
# gau_tol = [1.7e-3, 3e-4, 1e-5]
# gau_label = ["Gau. Loose", "Gau. Default", "Gau. Tight"]
# styles = ["-", "--", "-."]

# gau_tol = [1.7e-3, 3e-4]
# gau_label = ["Loose Tol.", " Default Tol.",]
# styles = ["-", "--"]


# for i, tol in enumerate(gau_tol):
#     plt.gca().axhline(
#         tol, label=gau_label[i], linewidth=2, linestyle=styles[i], color="k"
#     )

plt.xlim((epsilon1[0] * 1.1, epsilon1[-1]))
plt.xlabel(r"$\epsilon_1$ (Ha)")
plt.ylabel("Norm of Gradient Error (Ha/Bohr)")
plt.title(r"Gradient Error as a Function of $\epsilon_1$")
plt.ylim((1e-7, 3e-2))
plt.legend(framealpha=1.0)

plt.tight_layout()
plt.savefig("_figures/sc2_eps1_conv2.pdf")
plt.savefig("_figures/sc2_eps1_conv2.png", dpi=600)

