import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

df = pd.read_csv("uhf_survey.csv")

for obs in ["Energy (Ha)", "<S^2>"]:
    sns.set_context("talk")
    sns.catplot(
        x="Opt. Strategy",
        y=obs,
        row="Geometry",
        col="Multiplicity",
        hue="Stable",
        data=df,
        legend=True,
    )
    name = obs.split()[0]
    plt.savefig(f"_figures/uhf_survey_{name}.png")
    plt.close()
