import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.constants import physical_constants, N_A, calorie

ha_to_kcalmol = (
    physical_constants["hartree-joule relationship"][0] / calorie * N_A / 1e3
)

#
# Plot General Comparison
#

df = pd.read_csv("_data/dft_general.csv")


fig = plt.figure()
sns.set_style("ticks")
sns.set_palette(sns.color_palette(["#" + x for x in ["444444", "ff4343", "3725ff"]]))


# Plot energy
ax1 = plt.subplot(1, 2, 1)
df["Energy (Ha)"] = (
    np.array(
        [
            0,
            df.iloc[1]["Energy (Ha)"] - df.iloc[0]["Energy (Ha)"],
            df.iloc[2]["Energy (Ha)"] - df.iloc[0]["Energy (Ha)"],
            0,
            df.iloc[4]["Energy (Ha)"] - df.iloc[3]["Energy (Ha)"],
            df.iloc[5]["Energy (Ha)"] - df.iloc[3]["Energy (Ha)"],
            0,
            df.iloc[7]["Energy (Ha)"] - df.iloc[6]["Energy (Ha)"],
            df.iloc[8]["Energy (Ha)"] - df.iloc[6]["Energy (Ha)"],
        ]
    )
    * ha_to_kcalmol
)
sns.catplot(
    kind="point",
    data=df,
    x="Species",
    y="Energy (Ha)",
    hue="Functional",
    ax=ax1,
    legend=False,
)
ax1.set_ylabel("Relative Energy (kcal/mol)")
ax1.get_legend().remove()
plt.close()

# Plot Spin density
ax2 = plt.subplot(1, 2, 2)
sns.catplot(
    kind="point",
    data=df,
    x="Species",
    y="Fe Spin Density",
    hue="Functional",
    ax=ax2,
    # legend=False,
)
plt.close()


fig.tight_layout()
fig.savefig("_figures/dft_general_comparison.png", dpi=600)
fig.savefig("_figures/dft_general_comparison.pdf")

exit(0)

df = pd.read_csv("_data/dft_spin_data.csv")

plt.figure()
sns.set_style("ticks")
sns.set_palette(sns.color_palette(["#" + x for x in ["444444", "ff4343", "3725ff"]]))

sns.catplot(kind="bar", data=df, x="Species", y="Mulliken Spin Fe", hue="Functional")
plt.savefig("_figures/dft_spin_density.png", dpi=600)
plt.close()


sns.catplot(
    kind="bar", data=df, x="Species", y="$S^2$ Before Annhilation", hue="Functional"
)
plt.savefig("_figures/dft_spin_contamination.png", dpi=600)
plt.close()

#
# Plot geometry comparison
#


df = pd.read_csv("_data/dft_compare_funcs.csv")

plt.figure()
sns.set_style("ticks")
sns.set_palette(sns.color_palette(["#" + x for x in ["444444", "ff4343", "3725ff"]]))

sns.catplot(
    kind="bar",
    data=df,
    x="Species",
    y=r"Disp./atom from Ortu\~no (\AA)",
    hue="Functional",
)
plt.ylabel("Displacement/atom from Ortuño (Å)")

plt.savefig("_figures/dft_compare_geom.png", dpi=600)
plt.close()

df["Energy (Ha)"] = [
    0,
    df.iloc[1]["Energy (Ha)"] - df.iloc[0]["Energy (Ha)"],
    df.iloc[2]["Energy (Ha)"] - df.iloc[0]["Energy (Ha)"],
    0,
    df.iloc[4]["Energy (Ha)"] - df.iloc[3]["Energy (Ha)"],
    df.iloc[5]["Energy (Ha)"] - df.iloc[3]["Energy (Ha)"],
    0,
    df.iloc[7]["Energy (Ha)"] - df.iloc[6]["Energy (Ha)"],
    df.iloc[8]["Energy (Ha)"] - df.iloc[6]["Energy (Ha)"],
]
sns.catplot(
    kind="point", data=df, x="Species", y="Energy (Ha)", hue="Functional",
)

plt.savefig("_figures/dft_compare_energy.png", dpi=600)
plt.savefig("_figures/dft_compare_energy.pdf")
plt.close()

#
# Compare to Experimental
#
df = pd.read_csv("_data/dft_compare_funcs_geom.csv")
plt.figure()
sns.set_style("ticks")
sns.set_palette(sns.color_palette(["#" + x for x in ["011627", "f15152", "0044ff"]]))

df_bl = df.iloc[0:11]
df_ang = df.iloc[11:]
# df_bl.update(df.select_dtypes(include=[np.number]).abs())
# print(df_bl)
species = ["(1,A)*", "(3,A)*", "(3,C)*"] * 3
functionals = ["M06-2X"] * 3 + ["M06-L"] * 3 + ["MN15"] * 3
methods = [x + " " + y for x, y in zip(species, functionals)]
methods = methods[::3]

df_bl = pd.melt(
    df_bl,
    id_vars=["Structural Parameters"],
    value_vars=methods,
    var_name="Method",
    value_name=r"Error ($\AA$)",
)
# print(df_bl)
df_ang = pd.melt(
    df_ang,
    id_vars=["Structural Parameters"],
    value_vars=methods,
    var_name="Method",
    value_name="Error (Deg.)",
)

sns.catplot(
    data=df_bl,
    kind="bar",
    x="Structural Parameters",
    y=r"Error ($\AA$)",
    hue="Method",
    legend=False,
    aspect=1.5,
    # log=True,
)

for i in range(11):
    plt.vlines(i + 0.5, -1, 1, alpha=0.2)

plt.ylim((-0.05, 0.22))
plt.xticks(rotation=45, ha="right")
legend = plt.legend(loc="upper right")
plt.tight_layout()
plt.savefig("_figures/dft_compare_geom_bl.png", dpi=600)

