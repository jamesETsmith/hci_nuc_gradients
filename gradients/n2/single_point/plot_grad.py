import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#
# Read Data
#
casscf = np.load("_data/grad_casscf.npy")
vhciscf = np.load("_data/grad_vhciscf.npy")
vhciscf_aa = np.load("_data/grad_vhciscf_aa.npy")
hciscf = np.load("_data/grad_hciscf.npy")
hciscf_aa = np.load("_data/grad_hciscf_aa.npy")

n_atom = casscf.shape[0]
shci_grads = [vhciscf, vhciscf_aa, hciscf, hciscf_aa]
shci_errors = [np.linalg.norm(casscf - gr) for gr in shci_grads]

print(shci_errors)
print(casscf)
for gr in shci_grads:
    print(gr[0][2])


#
# Plot
#
tick_label = ["vHCISCF", "vHCISCF + AA", "HCISCF", "HCISCF + AA"]

plt.figure()
sns.set_style("darkgrid")
sns.set_palette("muted")
fs = 14  # Fontsize

sns.barplot(np.arange(len(shci_errors)), shci_errors)
plt.xticks(np.arange(len(tick_label)), tick_label)

plt.ylabel(r"$\mathcal{l}_1$ Norm of Gradient Error (Ha/Bohr)", fontsize=fs)
plt.title("Gradient Error as a Function of Method Variant", fontsize=fs * 1.2)
plt.tight_layout()
plt.savefig("_figures/n2_single_point_error.png", dpi=600)
plt.savefig("_figures/n2_single_point_error.pdf")

