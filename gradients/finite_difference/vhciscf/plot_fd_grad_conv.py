import json
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from pyscf.lib import param

# Collect data
hlist = [0.1, 0.01, 0.001, 0.0001]
data = {"h": [], "central_dif": [], "five_pt": [], "analytical": []}
frames = []
for i, h in enumerate(hlist):
    with open("_data/h={}.json".format(h), "r") as f:
        frames.append(pd.DataFrame(json.load(f), index=[i]))


df = pd.concat(frames)
df["error_central_dif"] = np.abs(df["central_dif"] - df["analytical"])
df["error_five_pt"] = np.abs(df["five_pt"] - df["analytical"])
df["h_bohr"] = df["h"] / param.BOHR
pd.set_option("precision", 8)
print(df)

# Plot
plt.figure()
sns.set_style("darkgrid")
plt.loglog(
    df["h_bohr"],
    df["error_central_dif"],
    marker="^",
    label=r"$|g_{central dif} - g_{analytical}|$",
)
plt.loglog(
    df["h_bohr"], df["error_five_pt"], marker="o", label=r"$|g_{5pt} - g_{analytical}|$"
)
plt.xlabel("Finite Difference Step Size (Bohr)")
plt.ylabel("Gradient Difference (Ha/Bohr)")
plt.title("vHCISCF")
plt.legend()
plt.savefig("_figures/fd_comparison.png", dpi=600)
