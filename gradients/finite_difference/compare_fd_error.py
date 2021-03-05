import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from pyscf.lib import param


# Collect data
hlist = [0.1, 0.01]
frames = []

for i, h in enumerate(hlist):
    with open("./casscf/_data/h={}.json".format(h), "r") as f:
        ddict = json.load(f)
        ddict["method"] = "CASSCF"
        frames.append(pd.DataFrame(ddict, index=[i]))

for i, h in enumerate(hlist):
    with open("./vhciscf/_data/h={}.json".format(h), "r") as f:
        ddict = json.load(f)
        ddict["method"] = "vHCISCF"
        frames.append(pd.DataFrame(ddict, index=[i]))

for i, h in enumerate(hlist):
    with open("./hciscf/_data/h={}.json".format(h), "r") as f:
        ddict = json.load(f)
        ddict["method"] = "HCISCF"
        frames.append(pd.DataFrame(ddict, index=[i]))


df = pd.concat(frames)
df["error_central_dif"] = np.abs(df["central_dif"] - df["analytical"])
df["error_five_pt"] = np.abs(df["five_pt"] - df["analytical"])
df["h_bohr"] = df["h"] / param.BOHR
pd.set_option("precision", 8)
print(df)


#
# Plot
#
sns.set_palette("muted")
sns.set_style("ticks")
sns.set_context("paper")

g = sns.catplot(
    x="h", y="error_five_pt", hue="method", data=df, kind="bar", palette="muted"
)
g._legend.remove()  # Remove Seaborn Legend so we can make our own
plt.gca().set(yscale="log")
plt.ylabel(r"| g$_{FD}$ - g$_{A}$ |")
plt.xlabel(r"Finite Difference Step ($\AA$)")
plt.title("Errors in Analytical Gradients\n(Compared to Finite Diff.)")
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)

plt.tight_layout()
plt.savefig("_figures/compare_fd_error.png", dpi=600)
plt.savefig("_figures/compare_fd_error.pdf")
