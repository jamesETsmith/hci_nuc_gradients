import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.constants import physical_constants, N_A, calorie

# Conversion factor
ha_to_kcalmol = (
    physical_constants["hartree-joule relationship"][0] / calorie * N_A / 1e3
)

os.makedirs("_data", exist_ok=True)

#
# Get Data
#
df = pd.read_csv("_data/40e_40o_initial_energy_comparison.csv")
df_1_A = pd.read_csv("_data/extrap_data/1_A.csv")
df_1_B = pd.read_csv("_data/extrap_data/1_B.csv")
df_1_C = pd.read_csv("_data/extrap_data/1_C.csv")
df_3_A = pd.read_csv("_data/extrap_data/3_A.csv")
df_3_B = pd.read_csv("_data/extrap_data/3_B.csv")
df_3_C = pd.read_csv("_data/extrap_data/3_C.csv")

#
# Print Summary Data
#
df_summary = pd.DataFrame(
    {
        "Multiplicity": [1, 1, 1, 3, 3, 3],
        "Geometry": ["A", "B", "C"] * 2,
        "Final Energy (Ha)": [
            df[(df["Geometry"] == "A") & (df["Multiplicity"] == 1)]["pt"].min(),
            df[(df["Geometry"] == "B") & (df["Multiplicity"] == 1)]["pt"].min(),
            df[(df["Geometry"] == "C") & (df["Multiplicity"] == 1)]["pt"].min(),
            df[(df["Geometry"] == "A") & (df["Multiplicity"] == 3)]["pt"].min(),
            df[(df["Geometry"] == "B") & (df["Multiplicity"] == 3)]["pt"].min(),
            df[(df["Geometry"] == "C") & (df["Multiplicity"] == 3)]["pt"].min(),
        ],
        "Extrap. Energy (Ha)": [
            df_1_A["Linear Fit Int. (Ha)"].min(),
            df_1_B["Linear Fit Int. (Ha)"].min(),
            df_1_C["Linear Fit Int. (Ha)"].min(),
            df_3_A["Linear Fit Int. (Ha)"].min(),
            df_3_B["Linear Fit Int. (Ha)"].min(),
            df_3_C["Linear Fit Int. (Ha)"].min(),
        ],
    }
)
print(df_summary)


#
# Plotting
#

plt.figure()
sns.set_style("ticks")
sns.set_palette("muted")
# sns.set_context("talk")

# sns.scatterplot(data=df, x="pt correction", y="pt", hue="Geometry", col="Multiplicity")
g = sns.FacetGrid(df, col="Multiplicity", hue="Geometry")
g.map(sns.scatterplot, "pt correction", "pt")
g.set_axis_labels("$E_2$ (Ha)", r"E$_{SHCI}$ (Ha)")

ax0 = g.facet_axis(0, 0)
ax1 = g.facet_axis(0, 1)


sns.lineplot(data=df_1_A, x="fit_x", y="lin_fit", ax=ax0)
sns.lineplot(data=df_1_B, x="fit_x", y="lin_fit", ax=ax0)
sns.lineplot(data=df_1_C, x="fit_x", y="lin_fit", ax=ax0)
sns.lineplot(data=df_3_A, x="fit_x", y="lin_fit", ax=ax1)
sns.lineplot(data=df_3_B, x="fit_x", y="lin_fit", ax=ax1)
sns.lineplot(data=df_3_C, x="fit_x", y="lin_fit", ax=ax1)

ax0.legend().set_visible(False)
ax1.legend().set_visible(False)

ax0.set_title("$S_z=0$")
ax1.set_title("$S_z=1$")

# plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
# plt.tight_layout()
g.add_legend()
# plt.tight_layout(pad=1)

os.makedirs("_figures", exist_ok=True)
plt.savefig(f"_figures/40e_40o_singet_energy_comparison.png", dpi=600)
plt.savefig(f"_figures/40e_40o_singet_energy_comparison.pdf")