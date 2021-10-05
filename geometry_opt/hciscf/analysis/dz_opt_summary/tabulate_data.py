import os
import numpy as np
import pandas as pd


def get_hciscf_data(logfile: str, noon_file: str, state: str) -> pd.Series:
    data = {
        "State": state,
        "<S^2>": None,
        "NO Energy (Ha)": None,
        "Energy (Ha)": None,
        "Fe Charge": None,
        "Fe 3dxy pop": None,
        "Fe 3dyz pop": None,
        "Fe 3dz^2 pop": None,
        "Fe 3dxz pop": None,
        "Fe 3dx2-y2 pop": None,
        "HONO-1": None,
        "HONO": None,
        "LUNO": None,
        "LUNO+1": None,
    }

    # Grab the NOONs
    if os.path.exists(noon_file):
        noons = np.loadtxt(noon_file)[117 : 117 + 4]
        data["HONO-1"] = noons[0]
        data["HONO"] = noons[1]
        data["LUNO"] = noons[2]
        data["LUNO+1"] = noons[3]
    else:
        data["HONO-1"] = None
        data["HONO"] = None
        data["LUNO"] = None
        data["LUNO+1"] = None

    # Parse log for data
    with open(logfile, "r") as f:
        lines = f.readlines()

    for line in lines:
        if "CASCI E" in line and "E(CI)" in line:
            lsplit = line.split()
            data["<S^2>"] = float(lsplit[9])
            data["Energy (Ha)"] = float(lsplit[3])
            # print(lsplit)
        if "CASCI energy" in line:
            lsplit = line.split()
            data["NO Energy (Ha)"] = float(lsplit[6])
            # print(lsplit)
        if "charge of  0" in line:
            lsplit = line.split()
            data["Fe Charge"] = float(lsplit[4])
        elif "pop of  0 Fe 3dxy" in line:
            data["Fe 3dxy pop"] = float(line.split()[-1])
        elif "pop of  0 Fe 3dyz" in line:
            data["Fe 3dyz pop"] = float(line.split()[-1])
        elif "pop of  0 Fe 3dz^2" in line:
            data["Fe 3dz^2 pop"] = float(line.split()[-1])
        elif "pop of  0 Fe 3dxz" in line:
            data["Fe 3dxz pop"] = float(line.split()[-1])
        elif "pop of  0 Fe 3dx2-y2" in line:
            data["Fe 3dx2-y2 pop"] = float(line.split()[-1])

        # Early stopping
        done = True
        for k, v in data.items():
            if v == None:
                done = False
                break

        if done:
            break

    return pd.Series(data)


#
# Main
#

os.makedirs("_initial_data", exist_ok=True)
os.makedirs("_final_data", exist_ok=True)

cas = "10e_10o"
# cas = "20e_20o"
# cas = "30e_30o"
cas = "40e_40o"
eps1 = "7.5e-05"
# eps1 = "0.0002"
states = [
    "1_A",
    "1_B",
    "1_C",
    "3_A",
    "3_B",
    "3_C",
]
noon_files = [
    f"../../dz/{d}/{cas}/_data/{d}_DZ_{cas}_eps1={eps1}_noons_initial.txt"
    for d in states
]
logfiles = [f"../../dz/{d}/{cas}/_logs/opt_{eps1}.out" for d in states]

df = pd.concat(
    [get_hciscf_data(*x) for x in zip(logfiles, noon_files, states)], axis=1
).T
print(df.to_string(float_format=(lambda x: f"{x:.6f}")))
df.to_csv(f"_initial_data/{cas}.csv", float_format=(lambda x: f"{x:.6f}"))

# df_noons = pd.DataFrame(np.array([np.loadtxt(nf) for nf in noon_files]).T, columns=states)
# df_noons.to_csv(f"_data/{cas}_noons.csv")

noon_files = [
    f"../../dz/{d}/{cas}/_data/{d}_DZ_{cas}_eps1={eps1}_noons_final.txt" for d in states
]
logfiles = [f"../../dz/{d}/{cas}/_logs/opt_{eps1}.out" for d in states]

df = pd.concat(
    [get_hciscf_data(*x) for x in zip(logfiles, noon_files, states)], axis=1
).T
print(df.to_string(float_format=(lambda x: f"{x:.6f}")))
df.to_csv(f"_final_data/{cas}.csv", float_format=(lambda x: f"{x:.6f}"))
