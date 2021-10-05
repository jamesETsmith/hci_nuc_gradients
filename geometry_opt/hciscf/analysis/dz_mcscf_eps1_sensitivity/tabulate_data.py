import sys
import os
import numpy as np
import pandas as pd

sys.path.append("..")
from analysis_utils import get_shci_data
from analysis_utils import get_natural_orbital_energy
from analysis_utils import get_hci_size

os.makedirs("_data", exist_ok=True)

# User settings
eps1s = [2e-4, 1e-4, 7.5e-5, 5e-5, 3e-5]
# eps1s = [1e-4, 7.5e-5, 5e-5]

# Grab energy data
mcscf_files = [f"../../dz/1_A/40e_40o_eps1_response/_logs/opt_{e}.out" for e in eps1s]
shci_files = [f"../../dz/1_A/40e_40o_eps1_response/_logs/{e}/output.dat" for e in eps1s]

mcscf_energies = [get_natural_orbital_energy(f) for f in mcscf_files]
mcscf_dets = [get_hci_size(f) for f in mcscf_files]
shci_energies = [get_shci_data(f)["pt"] for f in shci_files]
shci_error = [get_shci_data(f)["pt uncertainty"] for f in shci_files]
shci_dets = [get_hci_size(f) for f in shci_files]

# Gradients
#
eps1s = [2e-4, 1e-4, 7.5e-5, 5e-5, 3e-5]

gradients = [
    np.loadtxt(
        f"../../dz/1_A/40e_40o_eps1_response/_data/1_A_DZ_40e_40o_eps1={e}_gradients.txt"
    )
    for e in eps1s
]
n_atom_sqrt = np.sqrt(gradients[0].shape[0])
print(gradients[0].shape)
errors = [
    np.linalg.norm(gradients[i] - gradients[-1]) / n_atom_sqrt
    for i in range(len(gradients))
]
print(errors)
og = np.loadtxt("../../dz/1_A/40e_40o_eps1_response/_data/og.txt")
print(np.linalg.norm(og - gradients[-1]) / n_atom_sqrt)


# Setup and save dataframe
df = pd.DataFrame(
    {
        "Epsilon1 (Ha)": eps1s,
        "MCSCF Energy (Ha)": mcscf_energies,
        "MCSCF Determinants": mcscf_dets,
        "SHCI Energy (Ha)": shci_energies,
        "SHCI Energy Uncertainty (Ha)": shci_error,
        "SHCI Determinants": shci_dets,
        "Gradient Errors (Ha/Bohr)": errors,
    }
)

print(df)
df.to_csv("_data/mcscf_eps1_sensitivity.csv", index=False)

#
