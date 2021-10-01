import os
import numpy as np
import pandas as pd
import matplotlib
import seaborn as sns

plt = matplotlib.pyplot
os.makedirs("_figs", exist_ok=True)
matplotlib.rcParams["mathtext.fontset"] = "custom"
matplotlib.rcParams["mathtext.rm"] = "DejaVu Sans"
matplotlib.rcParams["mathtext.it"] = "DejaVu Sans"
matplotlib.rcParams["mathtext.bf"] = "DejaVu Sans:bold"
# sns.set_theme(context="paper", style="ticks", palette="muted")
sns.set_theme(style="ticks", palette="muted")


states = ["1_A", "1_B", "1_C", "3_A", "3_B", "3_C"]

# Plot the geometry errors for each state
for si in states:
    df = pd.read_csv(f"_data/{si}.csv")

    plt.figure()
    g = sns.barplot(x="Bond", y="Error", data=df, hue="Geometry")
    g.set_yscale("log")

    # Formatting
    plt.ylabel("Error $(\AA)$")
    plt.xticks(rotation=65, ha="right")
    plt.tight_layout()
    plt.savefig(f"_figs/{si}.png", dpi=600)
    plt.savefig(f"_figs/{si}.pdf")
    plt.close()
