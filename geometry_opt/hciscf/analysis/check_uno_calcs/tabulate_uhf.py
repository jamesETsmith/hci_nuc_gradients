import numpy as np
import pandas as pd


def get_scf_data(logfile: str, noon_file: str, state: str) -> pd.Series:
    data = {
        "State": state,
        "Stable": False,
        "Fe Spin Density": 0.0,
        "<S^2>": 0.0,
        "Energy (Ha)": 0.0,
        "HONO-1": 0.0,
        "HONO": 0.0,
        "LUNO": 0.0,
        "LUNO+1": 0.0,
    }

    # Parse data from logfile
    with open(logfile, "r") as f:
        lines = f.readlines()

    for line in lines:
        if "converged SCF" in line:
            lsplit = line.split()
            # print(lsplit)
            data["Energy (Ha)"] = float(lsplit[4])
            data["<S^2>"] = float(lsplit[7])
        if "spin density of  0" in line:
            lsplit = line.split()
            # print(lsplit)
            data["Fe Spin Density"] = float(lsplit[6])
        if "stability" in line:
            # print(line)
            if "wavefunction is stable" in line:
                data["Stable"] = True

    # Grab the frontier noons
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
noon_files = [f"../../{d}/gen_uno/_data/uhf_{d}_noons.txt" for d in states]
logfiles = [f"../../{d}/gen_uno/_logs/uhf_{d}.out" for d in states]

df = pd.concat([get_scf_data(*x) for x in zip(logfiles, noon_files, states)], axis=1,).T

print(df.to_string(float_format=(lambda x: f"{x:.5f}")))
df.to_csv("uhf.csv", float_format=(lambda x: f"{x:.6f}"))

