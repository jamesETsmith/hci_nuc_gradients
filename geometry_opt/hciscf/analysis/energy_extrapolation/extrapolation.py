import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from extrap_utils import extrapolate, get_shci_data
from scipy.constants import physical_constants, N_A, calorie

# Conversion factor
ha_to_kcalmol = (
    physical_constants["hartree-joule relationship"][0] / calorie * N_A / 1e3
)

os.makedirs("_data", exist_ok=True)

#
# CAS(40e,40o)
#

base_dir = "../../1_A/40e_40o/initial_extrap/"
i_1_A = extrapolate(base_dir, "1_A", "40e_40o", "initial")

base_dir = "../../1_A/40e_40o/final_extrap/"
f_1_A = extrapolate(base_dir, "1_A", "40e_40o", "final")

base_dir = "../../1_B/40e_40o/initial_extrap/"
i_1_B = extrapolate(base_dir, "1_B", "40e_40o", "initial")

base_dir = "../../1_B/40e_40o/final_extrap/"
f_1_B = extrapolate(base_dir, "1_B", "40e_40o", "final")

base_dir = "../../1_C/40e_40o/initial_extrap/"
i_1_C = extrapolate(base_dir, "1_C", "40e_40o", "initial")

base_dir = "../../1_C/40e_40o/final_extrap/"
f_1_C = extrapolate(base_dir, "1_C", "40e_40o", "final")

df_4040 = pd.DataFrame()
for d in [i_1_A, f_1_A, i_1_B, f_1_B, i_1_C, f_1_C]:
    df_4040 = df_4040.append(d, ignore_index=True)

df_4040 = df_4040[
    [
        "Multiplicity",
        "Geometry",
        "Initial/Final",
        "CAS",
        "Linear Fit Int. (Ha)",
        "Linear Fit Int. Uncertainty (Ha)",
        "Quadratic Fit Int. (Ha)",
        "Quadratic Fit Int. Uncertainty (Ha)",
    ]
]
print(df_4040)
df_4040.to_csv("_data/CAS(40e,40o)_extrap.csv", index=False)

exit(0)

#
# Extrapolated Gaps
#
# NOTE: Unit conversion to kcal/mol done at the very end
def calc_uncertainty(dE1: float, dE2: float) -> float:
    return np.sqrt(np.power(dE1, 2) + np.power(dE2, 2))


extr_gap = {
    "CAS": [
        "(40e,40o) Initial",
        "(40e,40o) Final",
        "(30e,30o) Initial TSG",
        "(30e,30o) Final TSG",
        "(40e,40o) Initial TSG",
        "(40e,40o) Final TSG",
    ],
    "Linear Fit Gap (kcal/mol)": [
        df_3030["Linear Fit Int. (Ha)"][2] - df_3030["Linear Fit Int. (Ha)"][0],
        df_3030["Linear Fit Int. (Ha)"][3] - df_3030["Linear Fit Int. (Ha)"][1],
        df_4040["Linear Fit Int. (Ha)"][2] - df_4040["Linear Fit Int. (Ha)"][0],
        df_4040["Linear Fit Int. (Ha)"][3] - df_4040["Linear Fit Int. (Ha)"][1],
        #
        df_3030["Linear Fit Int. (Ha)"][4] - df_3030["Linear Fit Int. (Ha)"][0],
        df_3030["Linear Fit Int. (Ha)"][5] - df_3030["Linear Fit Int. (Ha)"][1],
        df_4040["Linear Fit Int. (Ha)"][4] - df_4040["Linear Fit Int. (Ha)"][0],
        df_4040["Linear Fit Int. (Ha)"][5] - df_4040["Linear Fit Int. (Ha)"][1],
    ],
    "Linear Fit Gap Uncertainty (kcal/mol)": [
        calc_uncertainty(
            df_3030["Linear Fit Int. Uncertainty (Ha)"][0],
            df_3030["Linear Fit Int. Uncertainty (Ha)"][2],
        ),
        calc_uncertainty(
            df_3030["Linear Fit Int. Uncertainty (Ha)"][1],
            df_3030["Linear Fit Int. Uncertainty (Ha)"][3],
        ),
        calc_uncertainty(
            df_4040["Linear Fit Int. Uncertainty (Ha)"][0],
            df_4040["Linear Fit Int. Uncertainty (Ha)"][2],
        ),
        calc_uncertainty(
            df_4040["Linear Fit Int. Uncertainty (Ha)"][1],
            df_4040["Linear Fit Int. Uncertainty (Ha)"][3],
        ),
        #
        calc_uncertainty(
            df_3030["Linear Fit Int. Uncertainty (Ha)"][0],
            df_3030["Linear Fit Int. Uncertainty (Ha)"][4],
        ),
        calc_uncertainty(
            df_3030["Linear Fit Int. Uncertainty (Ha)"][1],
            df_3030["Linear Fit Int. Uncertainty (Ha)"][5],
        ),
        calc_uncertainty(
            df_4040["Linear Fit Int. Uncertainty (Ha)"][0],
            df_4040["Linear Fit Int. Uncertainty (Ha)"][4],
        ),
        calc_uncertainty(
            df_4040["Linear Fit Int. Uncertainty (Ha)"][1],
            df_4040["Linear Fit Int. Uncertainty (Ha)"][5],
        ),
    ],
    "Quadratic Fit Gap (kcal/mol)": [
        df_3030["Quadratic Fit Int. (Ha)"][2] - df_3030["Quadratic Fit Int. (Ha)"][0],
        df_3030["Quadratic Fit Int. (Ha)"][3] - df_3030["Quadratic Fit Int. (Ha)"][1],
        df_4040["Quadratic Fit Int. (Ha)"][2] - df_4040["Quadratic Fit Int. (Ha)"][0],
        df_4040["Quadratic Fit Int. (Ha)"][3] - df_4040["Quadratic Fit Int. (Ha)"][1],
        #
        df_3030["Quadratic Fit Int. (Ha)"][4] - df_3030["Quadratic Fit Int. (Ha)"][0],
        df_3030["Quadratic Fit Int. (Ha)"][5] - df_3030["Quadratic Fit Int. (Ha)"][1],
        df_4040["Quadratic Fit Int. (Ha)"][4] - df_4040["Quadratic Fit Int. (Ha)"][0],
        df_4040["Quadratic Fit Int. (Ha)"][5] - df_4040["Quadratic Fit Int. (Ha)"][1],
    ],
    "Quadratic Fit Gap Uncertainty (kcal/mol)": [
        calc_uncertainty(
            df_3030["Quadratic Fit Int. Uncertainty (Ha)"][0],
            df_3030["Quadratic Fit Int. Uncertainty (Ha)"][2],
        ),
        calc_uncertainty(
            df_3030["Quadratic Fit Int. Uncertainty (Ha)"][1],
            df_3030["Quadratic Fit Int. Uncertainty (Ha)"][3],
        ),
        calc_uncertainty(
            df_4040["Quadratic Fit Int. Uncertainty (Ha)"][0],
            df_4040["Quadratic Fit Int. Uncertainty (Ha)"][2],
        ),
        calc_uncertainty(
            df_4040["Quadratic Fit Int. Uncertainty (Ha)"][1],
            df_4040["Quadratic Fit Int. Uncertainty (Ha)"][3],
        ),
        #
        calc_uncertainty(
            df_3030["Quadratic Fit Int. Uncertainty (Ha)"][0],
            df_3030["Quadratic Fit Int. Uncertainty (Ha)"][4],
        ),
        calc_uncertainty(
            df_3030["Quadratic Fit Int. Uncertainty (Ha)"][1],
            df_3030["Quadratic Fit Int. Uncertainty (Ha)"][5],
        ),
        calc_uncertainty(
            df_4040["Quadratic Fit Int. Uncertainty (Ha)"][0],
            df_4040["Quadratic Fit Int. Uncertainty (Ha)"][4],
        ),
        calc_uncertainty(
            df_4040["Quadratic Fit Int. Uncertainty (Ha)"][1],
            df_4040["Quadratic Fit Int. Uncertainty (Ha)"][5],
        ),
    ],
}


df = pd.DataFrame(extr_gap)
df["Linear Fit Gap (kcal/mol)"] *= ha_to_kcalmol
df["Linear Fit Gap Uncertainty (kcal/mol)"] *= ha_to_kcalmol
df["Quadratic Fit Gap (kcal/mol)"] *= ha_to_kcalmol
df["Quadratic Fit Gap Uncertainty (kcal/mol)"] *= ha_to_kcalmol

print(df)
df.to_html("extrap_gaps.html", index=False)
df.to_csv("extrap_gaps.csv", index=False)
