import numpy as np
import pandas as pd


def get_hciscf_data(logfile: str, noon_file: str, state: str) -> pd.Series:
    data = {
        "State": state,
        # "Fe Spin Density": 0.0,
        "<S^2>": 0.0,
        "NO Energy (Ha)": 0.0,
        "Energy (Ha)": 0.0,
        "Fe Charge": 0.0,
        "Fe 3dxy pop": 0.0,
        "Fe 3dyz pop": 0.0,
        "Fe 3dz^2 pop": 0.0,
        "Fe 3dxz pop": 0.0,
        "Fe 3dx2-y2 pop": 0.0,
        "HONO-1": 0.0,
        "HONO": 0.0,
        "LUNO": 0.0,
        "LUNO+1": 0.0,
    }

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

    noons = np.loadtxt(noon_file)[117 : 117 + 4]
    data["HONO-1"] = noons[0]
    data["HONO"] = noons[1]
    data["LUNO"] = noons[2]
    data["LUNO+1"] = noons[3]

    return pd.Series(data)


states = [
    "1_A",
    "1_B",
    "1_C",
    "3_A",
    "3_B",
    "3_C",
]
noon_files = [
    f"../../tz/{d}/gen_uno/_data/{d}_TZ_100e_100o_eps1=0.0001_vhci_noons.txt"
    for d in states
]
logfiles = [
    f"../../tz/{d}/gen_uno/_logs/{d}_TZ_100e_100o_eps1=0.0001_vhci.out" for d in states
]

df = pd.concat(
    [get_hciscf_data(*x) for x in zip(logfiles, noon_files, states)], axis=1
).T
print(df.to_string(float_format=(lambda x: f"{x:.6f}")))
df.to_csv("uno_vhci.csv", float_format=(lambda x: f"{x:.6f}"))
