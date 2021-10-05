import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from extrap_utils import get_shci_data, extrapolate
from scipy.constants import physical_constants, N_A, calorie

# Conversion factor
ha_to_kcalmol = (
    physical_constants["hartree-joule relationship"][0] / calorie * N_A / 1e3
)

os.makedirs("_data", exist_ok=True)
os.makedirs("_data/extrap_data/", exist_ok=True)

#
# Get Data for
#

base_dir = "../../{}/40e_40o/initial_extrap/"
species = ["1_A", "1_B", "1_C", "3_A", "3_B", "3_C"]

df = pd.DataFrame()
for sp in species:
    df_species = pd.DataFrame()
    for i in range(10):
        path = os.path.join(base_dir.format(sp), f"output_{i}.dat")
        data = get_shci_data(path)
        data["Multiplicity"] = sp[0]
        data["Geometry"] = sp[2]
        data["pt correction"] = data["pt"] - data["var"]

        df = df.append(data, ignore_index=True)
        df_species = df_species.append(data, ignore_index=True)

    # Extrapolation
    extrap_data = extrapolate(df_species, sp, "40e_40o", "initial")
    pd.DataFrame(extrap_data).to_csv(f"_data/extrap_data/{sp}.csv", index=False)

print(df)

df.to_csv("_data/40e_40o_initial_energy_comparison.csv", index=False)