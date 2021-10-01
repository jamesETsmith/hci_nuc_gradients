import os
import sys
import numpy as np
import pandas as pd

sys.path.append("../")
from analysis_utils import parse_multi_XYZ

os.makedirs("_data", exist_ok=True)


def get_hciscf_data(xyz_file: str, noon_file: str, state: str) -> pd.Series:
    data = {
        "State": state,
        # "Fe Spin Density": 0.0,
        # "<S^2>": 0.0,
        # "NO Energy (Ha)": 0.0,
        "Energy (Ha)": 0.0,
        "Fe Charge": 0.0,
        "HONO-1": 0.0,
        "HONO": 0.0,
        "LUNO": 0.0,
        "LUNO+1": 0.0,
    }

    # Parse log for data
    energies, geometries, atoms = parse_multi_XYZ(xyz_file)
    data["Energy (Ha)"] = energies[0]

    noons = np.loadtxt(noon_file)[117 : 117 + 4]
    data["HONO-1"] = noons[0]
    data["HONO"] = noons[1]
    data["LUNO"] = noons[2]
    data["LUNO+1"] = noons[3]

    return pd.Series(data)


states = [
    "1_A",
    "1_B",
    "1_C",
    "3_A",
    "3_B",
    "3_C",
]


noon_files = [
    f"../../{d}/40e_40o/_data/{d}_40e_40o_eps1=0.0002_noons_initial.txt" for d in states
]
xyz_files = [f"../../{d}/40e_40o/_data/{d}_40e_40o_eps1=0.0002_optim.xyz" for d in states]
log_fils = [f""]

df = pd.concat(
    [get_hciscf_data(*x) for x in zip(xyz_files, noon_files, states)], axis=1
).T
print(df.to_string(float_format=(lambda x: f"{x:.6f}")))
df.to_csv("hciscf.csv", float_format=(lambda x: f"{x:.6f}"))
