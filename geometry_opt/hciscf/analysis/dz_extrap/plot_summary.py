import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append("../")
from plotting_utils import set_context, set_palette

df = pd.read_csv("_data/40e_40o_initial_extrap_comparison.csv")
df_full = pd.read_csv("_data/40e_40o_initial_energy_comparison_full.csv")

#
# Plot Intercept Comparison
#
plt.figure()
sns.set_palette("muted")
# sns.set_context("talk")

sns.scatterplot(
    x="Multiplicity",
    y="Linear Fit Int. (Ha)",
    hue="Geometry",
    # error="Linear Fit Int. Uncertainty (Ha)",
    data=df,
)

plt.ylabel("Energy (Ha)")
plt.xlabel("Mutliplicity")
plt.tight_layout()
plt.savefig("_figures/initial_summary.png", dpi=600)
plt.close()

#
# Plot Extrapolation Comparison
#
set_palette(3)
set_context("paper")
f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
# sns.set_palette("muted")
# sns.set_context("talk")


def linear_func(x, m, b):
    return x * m + b


def get_raw(df_full, m, g):
    df_subset = df_full.loc[df_full["Multiplicity"] == m]
    df_subset = df_subset.loc[df_subset["Geometry"] == g]
    # print(df_subset)
    x = df_subset["pt correction"]
    y = df_subset["pt"]
    return x, y


states = [(1, "A"), (1, "B"), (1, "C"), (3, "A"), (3, "B"), (3, "C")]


for i, state in enumerate(states):
    if i < 3:
        ax = ax1
    else:
        ax = ax2
    m, g = state
    x, y = get_raw(df_full, m, g)
    ax.scatter(x, y, label=f"{g}")

    # Fit
    fit_x = np.linspace(x.min(), 0)
    ax.plot(
        fit_x, linear_func(fit_x, df["Linear Slope"][i], df["Linear Fit Int. (Ha)"][i])
    )


ax1.set_xlabel("$E_2$ (Ha)")
ax2.set_xlabel("$E_2$ (Ha)")
ax1.set_ylabel(r"E$_{SHCI}$ (Ha)")
ax1.set_xlim(-0.03, 0.001)
ax2.set_xlim(-0.03, 0.001)
ax1.get_yaxis().get_major_formatter().set_useOffset(False)
ax1.get_yaxis().get_major_formatter().set_scientific(False)
ax1.legend()

ax1.set_title("Singlet")
ax2.set_title("Triplet")

plt.tight_layout()
plt.savefig("_figures/initial_extrapolation_summary.png", dpi=600)
plt.savefig("_figures/initial_extrapolation_summary.pdf")
plt.close()

#
# Plot Final Comparison
#
df = pd.read_csv("_data/40e_40o_final_extrap_comparison.csv")
df_full = pd.read_csv("_data/40e_40o_final_energy_comparison_full.csv")

set_palette(3)
set_context("paper")
f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
# sns.set_palette("muted")
# sns.set_context("talk")


def linear_func(x, m, b):
    return x * m + b


def get_raw(df_full, m, g):
    df_subset = df_full.loc[df_full["Multiplicity"] == m]
    df_subset = df_subset.loc[df_subset["Geometry"] == g]
    # print(df_subset)
    x = df_subset["pt correction"]
    y = df_subset["pt"]
    return x, y


states = [(1, "A"), (1, "B"), (1, "C"), (3, "A"), (3, "B"), (3, "C")]


for i, state in enumerate(states):
    if state[0] == 1:
        ax = ax1
    else:
        ax = ax2
    m, g = state
    x, y = get_raw(df_full, m, g)
    ax.scatter(x, y, label=f"{g}")

    # Fit
    fit_x = np.linspace(x.min(), 0)
    ax.plot(
        fit_x, linear_func(fit_x, df["Linear Slope"][i], df["Linear Fit Int. (Ha)"][i])
    )


ax1.set_xlabel("$E_2$ (Ha)")
ax2.set_xlabel("$E_2$ (Ha)")
ax1.set_ylabel(r"E$_{SHCI}$ (Ha)")
ax1.set_xlim(-0.03, 0.001)
ax2.set_xlim(-0.03, 0.001)
ax1.get_yaxis().get_major_formatter().set_useOffset(False)
ax1.get_yaxis().get_major_formatter().set_scientific(False)
ax1.legend()

ax1.set_title("Singlet")
ax2.set_title("Triplet")

plt.tight_layout()
plt.savefig("_figures/final_extrapolation_summary.png", dpi=600)
plt.savefig("_figures/final_extrapolation_summary.pdf")
plt.close()
