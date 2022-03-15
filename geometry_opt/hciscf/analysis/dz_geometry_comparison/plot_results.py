import os, sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sys.path.append("../")
from plotting_utils import set_context, set_palette

os.makedirs("_figs", exist_ok=True)


set_context("paper")
set_palette(4)

states = ["1_A", "1_B", "1_C", "3_A", "3_B", "3_C"]

# Plot the geometry errors for each state
for si in states:
    df = pd.read_csv(f"_data/{si}.csv")

    plt.figure()
    g = sns.barplot(x="Bond", y="Error", data=df, hue="Method")
    g.set_yscale("log")

    # Formatting
    plt.ylabel("Error $(\AA)$")
    plt.xticks(rotation=65, ha="right")
    plt.tight_layout()
    plt.savefig(f"_figs/{si}.png", dpi=600)
    plt.savefig(f"_figs/{si}.pdf")
    plt.close()
