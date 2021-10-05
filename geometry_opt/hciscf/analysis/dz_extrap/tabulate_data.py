import os, sys
import numpy as np
import pandas as pd
from scipy.constants import physical_constants, N_A, calorie

sys.path.append("..")
from analysis_utils import get_shci_data, extrapolate

# Conversion factor
ha_to_kcalmol = (
    physical_constants["hartree-joule relationship"][0] / calorie * N_A / 1e3
)

os.makedirs("_data", exist_ok=True)
os.makedirs("_data/extrap_data/", exist_ok=True)

#
# Get Data for initial
#

base_dir = "../../dz/{}/40e_40o/initial_extrapolation/"
species = ["1_A", "1_B", "1_C", "3_A", "3_B", "3_C"]
# species = ["1_A", "3_A", "3_B"]

df_full = pd.DataFrame()
df_extrap = pd.DataFrame()
for sp in species:
    df_species = pd.DataFrame()
    for i in range(10):
        path = os.path.join(base_dir.format(sp), f"output_{i}.dat")
        data = get_shci_data(path)
        data["Multiplicity"] = sp[0]
        data["Geometry"] = sp[2]
        data["pt correction"] = data["pt"] - data["var"]

        df_full = df_full.append(data, ignore_index=True)
        df_species = df_species.append(data, ignore_index=True)

    # Extrapolation
    print(df_species)
    if not np.any(df_species["pt uncertainty"] == 1.0):
        extrap_data = extrapolate(df_species, sp, "40e_40o", "initial")
        print(extrap_data)
        pd.Series(extrap_data).to_csv(f"_data/extrap_data/{sp}.csv")
        df_extrap = df_extrap.append(pd.Series(extrap_data), ignore_index=True)

# print(df)

df_full.to_csv("_data/40e_40o_initial_energy_comparison_full.csv", index=False)
df_extrap.to_csv("_data/40e_40o_initial_extrap_comparison.csv", index=False)

#
# Get Data for final
#

base_dir = "../../dz/{}/40e_40o/final_extrapolation/"
species = ["1_A", "1_B", "1_C", "3_A", "3_B", "3_C"]
# species = ["1_A", "3_A", "3_B"]

df_full = pd.DataFrame()
df_extrap = pd.DataFrame()
for sp in species:
    df_species = pd.DataFrame()
    for i in range(10):
        path = os.path.join(base_dir.format(sp), f"output_{i}.dat")
        data = get_shci_data(path)
        data["Multiplicity"] = sp[0]
        data["Geometry"] = sp[2]
        data["pt correction"] = data["pt"] - data["var"]

        df_full = df_full.append(data, ignore_index=True)
        df_species = df_species.append(data, ignore_index=True)

    # Extrapolation
    print(df_species)
    if not np.any(df_species["pt uncertainty"] == 1.0):
        extrap_data = extrapolate(df_species, sp, "40e_40o", "final")
        print(extrap_data)
        pd.Series(extrap_data).to_csv(f"_data/extrap_data/{sp}.csv")
        df_extrap = df_extrap.append(pd.Series(extrap_data), ignore_index=True)

# print(df)

df_full.to_csv("_data/40e_40o_final_energy_comparison_full.csv", index=False)
df_extrap.to_csv("_data/40e_40o_final_extrap_comparison.csv", index=False)
