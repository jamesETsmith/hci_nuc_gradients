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
sns.set_theme(context="talk", style="ticks", palette="muted")

# Plot Energies
plt.figure()

# Plot actual data
plt.plot(df["Epsilon1 (Ha)"], df["Gradient Errors (Ha/Bohr)"], "o-")

# Plot chemical accuracy
best_answer = df["Gradient Errors (Ha/Bohr)"].iloc[-1]
tol = 3e-4  # Ha/Bohr
upper_bound = np.ones(df["Epsilon1 (Ha)"].size) * best_answer + tol
lower_bound = np.ones(df["Epsilon1 (Ha)"].size) * best_answer
plt.gca().fill_between(
    df["Epsilon1 (Ha)"],
    lower_bound,
    upper_bound,
    facecolor="gray",
    alpha=0.5,
    label="Tolerance",
)

# Formatting
plt.legend()
plt.xticks(df["Epsilon1 (Ha)"][::2], [f"{e:.1e}" for e in df["Epsilon1 (Ha)"][::2]])
# plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
# plt.gca().get_yaxis().get_major_formatter().set_scientific(False)
plt.xlabel("$\epsilon_1 (Ha)$")
plt.ylabel("Gradient Error (Ha/Bohr)")
plt.tight_layout()
plt.savefig("_figs/grad_vs_eps1.png", dpi=600)
plt.savefig("_figs/grad_vs_eps1.pdf")
plt.close()
