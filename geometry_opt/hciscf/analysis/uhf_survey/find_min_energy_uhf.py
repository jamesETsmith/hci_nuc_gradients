import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Gather data
df_diis = pd.read_csv("_data/uhf_survey_dft_no_adiis.csv", index_col=0)
df_adiis = pd.read_csv("_data/uhf_survey.csv", index_col=0)
df_newton = pd.read_csv("_data/uhf_survey_dft_newton.csv", index_col=0)
df_nofn_diis = pd.read_csv("_data/uhf_survey_no_fn_dft_diis.csv", index_col=0)
df_nofn_adiis = pd.read_csv("_data/uhf_survey_no_fn_dft_adiis.csv", index_col=0)


# Add columns for DFT Opt Strategy and combine dataframes
dfs = [df_diis, df_adiis, df_nofn_diis, df_nofn_adiis, df_newton]
dft_opt_strategy = ["Fast_Newton+DIIS", "Fast_Newton+ADIIS", "DIIS", "ADIIS", "NEWTON"]

for i, df in enumerate(dfs):
    df["DFT Opt. Strategy"] = [dft_opt_strategy[i]] * len(df.index)

df = pd.concat(dfs, ignore_index=True)
# Cleanup
df["Opt. Strategy"] = [x.upper() for x in df["Opt. Strategy"]]
df = df.rename(columns={"Opt. Strategy": "UHF Opt. Strategy"})
print(df)

# Search combined dataframe for lowest energy solution for each state
multiplicity = [1, 3]
geometry = ["A", "B", "C"]

df_most_stable = pd.DataFrame()
os.makedirs("_data/by_species/", exist_ok=True)

for m in multiplicity:
    for g in geometry:
        df_m_g = df[(df["Multiplicity"] == m) & (df["Geometry"] == g)].reset_index(
            drop=True
        )
        if len(df_m_g) != 15:
            raise AssertionError(
                f"Something went wrong with selecting data for species {m}_{g}"
            )
        print(df_m_g.sort_values(by="Energy (Ha)"))
        most_stable = df_m_g.iloc[df_m_g["Energy (Ha)"].idxmin()]
        df_m_g.sort_values(by="Energy (Ha)").to_csv(f"_data/by_species/{m}_{g}.csv")
        print(most_stable)
        df_most_stable = pd.concat(
            [df_most_stable, most_stable], axis=1, ignore_index=True
        )

df_most_stable = df_most_stable.T
print(df_most_stable)
df_most_stable.to_csv("_data/uhf_survey_most_stable.csv")

#
sns.set_context("talk")
sns.set_style("ticks")
g = sns.catplot(
    data=df,
    x="UHF Opt. Strategy",
    y="Energy (Ha)",
    row="Multiplicity",
    col="Geometry",
    hue="DFT Opt. Strategy",
)
plt.ticklabel_format(axis="y", useOffset=False)
g.fig.subplots_adjust(left=0.12)
plt.savefig("_figures/uhf_most_stable.png", dpi=600)
