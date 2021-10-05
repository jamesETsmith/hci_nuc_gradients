import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def calc_ncas(noons: np.ndarray, tol: float) -> int:
    cas_idx = np.where((noons > tol) & (noons < (2.0 - tol)))
    # print(cas_idx)
    return cas_idx[0].size


df = pd.read_csv("noons.csv")


states = [
    "1_A",
    "1_B",
    "1_C",
    "3_A",
    "3_B",
    "3_C",
]
#
# Plot
#
plt.figure()
for s in states:
    noons = df[s]
    # noons = np.array([min(abs(noon-1.5), abs(noon-0.5)) for noon in noons])
    plt.plot(np.arange(noons.size) + 1e-14, noons, label=s)

plt.xlim((90, 140))
plt.legend()
plt.savefig("noons.png")
plt.close()

# Plot
tols = np.logspace(-1, -2.5, num=50)
print(tols)
plt.figure(figsize=(12, 6))
data = {}
for s in states:
    ncas_sizes = [calc_ncas(df[s], t) for t in tols]
    print(ncas_sizes)
    plt.semilogx(tols, ncas_sizes, "o-", label=s)
plt.legend()
plt.savefig("cas_sizes.png")
plt.close()

# Plot for CCQ retreat
plt.figure()
noons = df["1_A"]
# noons = np.array([min(abs(noon - 1.5), abs(noon - 0.5)) for noon in noons])
plt.plot(np.arange(noons.size), noons, "o-")

plt.xlim((90, 140))
plt.ylim((1, 2))
plt.savefig("noons_ccq_retreat.png", dpi=600)
plt.savefig("noons_ccq_retreat.pdf")
plt.close()
