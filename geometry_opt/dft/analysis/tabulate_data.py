"""
Collect the following data from the DFT results:

- Energies
- Geometries
- Natural orbital occupation numbers (NOONs)
- Spin densities on Fe
- Spin contamination


Compare a subset of the strucutral parameters to those reported by Stieber.
 
"""

import os
import numpy as np
import pandas as pd
import cclib
from scipy.constants import physical_constants, N_A, calorie
from rmsd import kabsch_rmsd

# Helper functions
def format_exponential(x: float) -> str:
    """Format exponential float values to look pretty in LaTeX.
    """
    if x == 0:
        return "0.0"

    float_str = f"{x:e}"
    base, exponent = float_str.replace("e-0", "e-").replace("e+0", "e+").split("e")
    return r"${0:.1f} \cdot 10^{{{1}}}$".format(float(base), exponent)


def get_ss_and_s(filename: str) -> dict:

    with open(filename, "r") as f:
        lines = f.readlines()

    for l in reversed(lines):
        if "<S**2" in l:
            lsplit = l.split()
            print(lsplit)
            ss = float(lsplit[7])
            s = float(lsplit[9])
            return ss, s


# Conversion factors
ha_to_ev = physical_constants["Hartree energy in eV"][0]
ha_to_kcalmol = (
    physical_constants["hartree-joule relationship"][0] / calorie * N_A / 1e3
)


#
# Set up output directories
#
for d in ["_data", "_figures", "_tables"]:
    os.makedirs(d, exist_ok=True)

#
# Collect General Data
#

g16_dirs = [
    "M062X_1_A",
    "M062X_3_A",
    "M062X_3_C",
    "M06L_1_A",
    "M06L_3_A",
    "M06L_3_C",
    "MN15_1_A",
    "MN15_3_A",
    "MN15_3_C",
]


energies = []
geometries = []
# species = []
multiplicity = []
starting_geom = []
functional = []
fe_spin = []
spin = []
spin_squared = []
geom_11_bs = None  # Using the notation of Ortuno and Cramer
geom_13_bs = None  # Using the notation of Ortuno and Cramer

# Get energies, species, and geometries
for d in g16_dirs:
    path = os.path.join("..", d, "_opt.log")
    ccdata = cclib.io.ccread(path)

    # Get functional and species info
    func = "M06-L"
    if "M062X" in d:
        func = "M06-2X"
    elif "MN15" in d:
        func = "MN15"
    functional.append(func)
    # species.append(f"$({d[-3]},{d[-1]})^*$")
    multiplicity.append(d[-3])
    starting_geom.append(r"\textbf{" + d[-1] + r"}")

    # Get general info
    geometries.append(ccdata.atomcoords[-1])
    energies.append(ccdata.scfenergies[-1] / ha_to_ev)
    fe_spin.append(abs(ccdata.atomspins["mulliken"][0]))
    ss, s = get_ss_and_s(path)
    spin.append(s)
    spin_squared.append(ss)

    if d == "M06L_1_A":
        geom_11_bs = ccdata.atomcoords[0]
    elif d == "M06L_3_C":
        geom_13_bs = ccdata.atomcoords[0]


# Get RMSD compared to Ortuno and Cramer
displacement = []

for i, gi in enumerate(geometries):
    comp_geom = geom_11_bs
    # Compare triplet geometries to Ortuno and Cramer's triplet
    if i % 3 == 1 or i % 3 == 2:
        comp_geom = geom_13_bs
    displacement.append(kabsch_rmsd(gi, comp_geom, translate=True))

# Collect data into a Pandas DataFrame
df_general = pd.DataFrame(
    {
        "Functional": functional,
        # "Species": species,
        "Multiplicity": multiplicity,
        "Starting Geometry": starting_geom,
        "Energy (Ha)": energies,
        "Fe Spin Density": fe_spin,
        # r"$\ev{\hat{S}}$": spin,
        r"$\ev{\hat{S}^2}$": spin_squared,
        "RMSD (Å)": displacement,
    }
)

df_geom = pd.DataFrame(
    {
        "Functional": functional,
        "Multiplicity": multiplicity,
        "Starting Geometry": starting_geom,
        "Geometry (Å)": geometries,
    }
)

# Save dataframes
df_general.to_csv("_data/dft_general.csv")
df_geom.to_csv("_data/dft_geoms.csv")

# Save tables for latex
formatters = {
    "Energy (Ha)": lambda x: f"{x:.6f}",
    "RMSD (Å)": format_exponential,
}
df_general.round(12).to_latex(
    buf="_tables/dft_general_comp.tex",
    index=False,
    escape=False,
    formatters=formatters,
    column_format="c" * len(df_general.columns),
)


#
# Collect NOONs
#

noons = []
for d in g16_dirs:
    path = os.path.join("..", d, "_data", f"{d}_noons.txt")
    noons.append(np.loadtxt(path))

noons = np.array(noons)

df_noons = pd.DataFrame(
    noons.T,
    columns=[
        x + " " + y + " " + z
        for x, y, z in zip(functional, multiplicity, starting_geom)
    ],
)
df_noons.to_csv("_data/dft_noons.csv")

frontier_noons = noons[:, 117:121]
df_frontier_noons = pd.DataFrame(
    frontier_noons, columns=["HONO-1", "HONO", "LUNO", "LUNO+1"]
)
df_frontier_noons["Functional"] = functional
# df_frontier_noons["Species"] = species
df_frontier_noons["Multiplicity"] = multiplicity
df_frontier_noons["Starting Geometry"] = starting_geom

df_frontier_noons = df_frontier_noons.reindex(
    [
        "Functional",
        "Multiplicity",
        "Starting Geometry",
        "HONO-1",
        "HONO",
        "LUNO",
        "LUNO+1",
    ],
    axis=1,
)

df_frontier_noons.to_csv("_data/dft_frontier_noons.csv")
print(df_frontier_noons)
df_frontier_noons.round(12).to_latex(
    buf="_tables/dft_frontier_noons.tex",
    index=False,
    escape=False,
    float_format=lambda x: f"{x:.3f}",
    column_format="c" * len(df_frontier_noons.columns),
)

print("STOPPING EARLY")
exit(0)

#
# Stieber Data
#
stieber_data = [
    # In angstroms
    ["Fe(1)−N(1)", 1.914, 0, 2],
    ["Fe(1)−N(2)", 1.842, 0, 1],
    ["Fe(1)−N(3)", 1.904, 0, 3],
    ["Fe(1)−N(4)", 1.799, 0, 4],
    ["N(4)−N(5)", 1.117, 4, 5],
    ["N(1)−C(2)", 1.344, 1, 6],
    ["N(3)−C(8)", 1.342, 3, 9],
    ["N(2)−C(3)", 1.370, 2, 7],
    ["N(2)−C(7)", 1.372, 2, 8],
    ["C(2)−C(3)", 1.440, 6, 7],
    ["C(7)−C(8)", 1.432, 8, 9],
    # In Degrees
    ["N(1)−Fe(1)−N(2)", 80.82, 1, 0, 2],
    ["N(1)−Fe(1)−N(3)", 160.71, 1, 0, 3],
    ["N(1)−Fe(1)−N(4)", 98.81, 1, 0, 4],
    ["N(2)−Fe(1)−N(3)", 80.64, 2, 0, 3],
    ["N(2)−Fe(1)−N(4)", 169.42, 2, 0, 4],
    ["N(3)−Fe(1)−N(4)", 98.44, 3, 0, 4],
]

#
# Helpers
#
def calculate_bl(indices: list, geom: np.ndarray) -> float:
    atom1 = geom[indices[0], :]
    atom2 = geom[indices[1], :]
    return np.linalg.norm(np.array(atom1) - np.array(atom2))


def calculate_bond_angle(indices: list, geom: np.ndarray) -> float:
    atom1 = geom[indices[0], :]
    center = geom[indices[1], :]
    atom2 = geom[indices[2], :]
    v1 = np.array(atom1) - np.array(center)
    v2 = np.array(atom2) - np.array(center)
    alpha = np.arccos(v1.dot(v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    return 180 * alpha / np.pi


#
# Compare our geometry data to Stieber's reported structural parameters
#


data = {"Method": ["Expt.", "Ortuno"]}
for i in stieber_data:
    data[i[0]] = [i[1]]

    # Bond length
    if len(i) == 4:
        data[i[0]].append(calculate_bl(i[2:], geom_11_bs))
    elif len(i) == 5:
        data[i[0]].append(calculate_bond_angle(i[2:], geom_11_bs))


for i, g in enumerate(geometries):
    data["Method"].append(species[i])
    for sd in stieber_data:
        # Bond length
        if len(sd) == 4:
            data[sd[0]].append(calculate_bl(sd[2:], g))
        elif len(sd) == 5:
            data[sd[0]].append(calculate_bond_angle(sd[2:], g))

df = pd.DataFrame(data).set_index("Method").T
df.index.name = "Structural Parameters"
print(df)
df = df.sub(df["Expt."], axis="index")
print(df.to_string(float_format=(lambda x: f"{x:.1e}")))
df.to_csv("_data/dft_compare_to_expt.csv")

