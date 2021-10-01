import os
import numpy as np
import pandas as pd
import matplotlib
import seaborn as sns

# The the same font throughout
plt = matplotlib.pyplot
matplotlib.rcParams["mathtext.fontset"] = "custom"
matplotlib.rcParams["mathtext.rm"] = "DejaVu Sans"
matplotlib.rcParams["mathtext.it"] = "DejaVu Sans"
matplotlib.rcParams["mathtext.bf"] = "DejaVu Sans:bold"

# Setup output dirs
os.makedirs("_figs", exist_ok=True)

# Read in data
df = pd.read_csv("_data/mcscf_eps1_sensitivity.csv")
print(df)

# Set theme
# sns.set_theme(context="talk", style="ticks", palette="muted")
sns.set_theme(style="ticks", palette="muted")

# Plot Energies
plt.figure()

# Plot actual data
plt.plot(df["Epsilon1 (Ha)"], df["MCSCF Energy (Ha)"], "o-", label="MCSCF")
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
plt.semilogy(df["Epsilon1 (Ha)"], df["MCSCF Determinants"], "o-", label="MCSCF")
plt.semilogy(df["Epsilon1 (Ha)"], df["SHCI Determinants"], "o-", label="SHCI")

# Formatting
plt.legend()
plt.xticks(df["Epsilon1 (Ha)"][::2], [f"{e:.1e}" for e in df["Epsilon1 (Ha)"][::2]])
plt.ylim((10 ** 6, 10 ** 8))
plt.xlabel("$\epsilon_1 (Ha)$")
plt.ylabel("# Determinants")
plt.tight_layout()
plt.savefig("_figs/mcscf_dets_vs_eps1.png", dpi=600)
plt.savefig("_figs/mcscf_dets_vs_eps1.pdf")
plt.close()
