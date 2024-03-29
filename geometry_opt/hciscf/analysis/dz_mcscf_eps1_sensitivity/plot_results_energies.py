import os
import numpy as np
import pandas as pd
import matplotlib
import seaborn as sns
import matplotlib.pyplot as plt

import sys

sys.path.append("..")
from plotting_utils import set_context, set_palette

# The the same font throughout


# Setup output dirs
os.makedirs("_figs", exist_ok=True)

# Read in data
df = pd.read_csv("_data/mcscf_eps1_sensitivity.csv")
print(df)

# Set theme
set_palette(2)
set_context("paper", font_scale=1.25)

# Plot Energies
plt.figure()

# Plot actual data
plt.plot(df["Epsilon1 (Ha)"], df["MCSCF Energy (Ha)"], "o-", label="vHCISCF")
plt.plot(df["Epsilon1 (Ha)"], df["SHCI Energy (Ha)"], "o-", label="SHCI")

# Plot chemical accuracy
best_answer = df["SHCI Energy (Ha)"].iloc[-1]
chemical_accuracy = 0.0015936  # Ha
upper_bound = np.ones(df["Epsilon1 (Ha)"].size) * best_answer + chemical_accuracy
lower_bound = np.ones(df["Epsilon1 (Ha)"].size) * best_answer - chemical_accuracy
plt.gca().fill_between(
    df["Epsilon1 (Ha)"],
    lower_bound,
    upper_bound,
    facecolor="gray",
    alpha=0.5,
    label="Chemical Acc.",
)

# Formatting
plt.legend()
plt.xticks(df["Epsilon1 (Ha)"][::2], [f"{e:.1e}" for e in df["Epsilon1 (Ha)"][::2]])
plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
plt.gca().get_yaxis().get_major_formatter().set_scientific(False)
plt.xlabel("$\epsilon_1 (Ha)$")
plt.ylabel("Energy (Ha)")
plt.tight_layout()
plt.savefig("_figs/mcscf_eps1_sensitivity.png", dpi=600)
plt.savefig("_figs/mcscf_eps1_sensitivity.pdf")
plt.close()

# Plot determinants
plt.figure()

# Plot actual data
plt.semilogy(df["Epsilon1 (Ha)"], df["MCSCF Determinants"], "o-", label="vHCISCF")
plt.semilogy(df["Epsilon1 (Ha)"], df["SHCI Determinants"], "o-", label="SHCI")

# Formatting
plt.legend()
plt.xticks(df["Epsilon1 (Ha)"][::2], [f"{e:.1e}" for e in df["Epsilon1 (Ha)"][::2]])
plt.ylim((10 ** 6, 10 ** 8))
plt.xlabel("$\epsilon_1 (Ha)$")
plt.ylabel("Number of Determinants")
plt.tight_layout()
plt.savefig("_figs/mcscf_dets_vs_eps1.png", dpi=600)
plt.savefig("_figs/mcscf_dets_vs_eps1.pdf")
plt.close()
