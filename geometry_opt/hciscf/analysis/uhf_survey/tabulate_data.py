import os
import numpy as np
import pandas as pd

#
# Set up dirs
#
for d in ["_data", "_figures"]:
    os.makedirs(d, exist_ok=True)


def get_scf_data(
    logfile: str, noon_file: str, state: str, opt_strategy: str
) -> pd.Series:
    data = {
        "Multiplicity": state[0],
        "Geometry": state[-1],
        "Opt. Strategy": opt_strategy,
        "SCF Converged": False,
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

    for i, line in enumerate(lines):
        if "converged SCF" in line:
            lsplit = line.split()
            # print(lsplit)
            data["Energy (Ha)"] = float(lsplit[4])
            data["<S^2>"] = float(lsplit[7])
            data["SCF Converged"] = True
        if "SCF not converged." in line:
            lsplit = lines[i + 1].split()
            data["Energy (Ha)"] = float(lsplit[3])
            data["<S^2>"] = float(lsplit[9])
            data["SCF Converged"] = False

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
opt_stategies = ["adiis", "diis", "newton"]


#
# Get data for original survey
#

series = []
for s in states:
    for opt in opt_stategies:
        noon_file = f"../../uhf_survey/{s}/_data/uhf_{s}_{opt}_noons.txt"
        log_file = f"../../uhf_survey/{s}/_logs/uhf_{s}_{opt}.out"
        series.append(get_scf_data(log_file, noon_file, s, opt))


df = pd.concat(series, axis=1,).T
print(df)
df.to_csv("_data/uhf_survey.csv")

#
# Get data for survey where DFT calc uses CDIIS (rather than ADIIS)
#
series = []
for s in states:
    for opt in opt_stategies:
        noon_file = f"../../uhf_survey_dft_no_adiis/{s}/_data/uhf_{s}_{opt}_noons.txt"
        log_file = f"../../uhf_survey_dft_no_adiis/{s}/_logs/uhf_{s}_{opt}.out"
        series.append(get_scf_data(log_file, noon_file, s, opt))


df = pd.concat(series, axis=1,).T
print(df)
df.to_csv("_data/uhf_survey_dft_no_adiis.csv")


#
# Get data for survey where DFT calc uses Newton (rather than ADIIS)
#
series = []
for s in states:
    for opt in opt_stategies:
        noon_file = f"../../uhf_survey_dft_newton/{s}/_data/uhf_{s}_{opt}_noons.txt"
        log_file = f"../../uhf_survey_dft_newton/{s}/_logs/uhf_{s}_{opt}.out"
        series.append(get_scf_data(log_file, noon_file, s, opt))


df = pd.concat(series, axis=1,).T
print(df)
df.to_csv("_data/uhf_survey_dft_newton.csv")


#
# Get data for survey where DFT calc uses ADIIS alone rather than fast_newton
#
series = []
for s in states:
    for opt in opt_stategies:
        noon_file = (
            f"../../uhf_survey_no_fn_dft_adiis/{s}/_data/uhf_{s}_{opt}_noons.txt"
        )
        log_file = f"../../uhf_survey_no_fn_dft_adiis/{s}/_logs/uhf_{s}_{opt}.out"
        series.append(get_scf_data(log_file, noon_file, s, opt))


df = pd.concat(series, axis=1,).T
print(df)
df.to_csv("_data/uhf_survey_no_fn_dft_adiis.csv")

#
# Get data for survey where DFT calc uses DIIS w/o fast newton
series = []
for s in states:
    for opt in opt_stategies:
        noon_file = f"../../uhf_survey_no_fn_dft_diis/{s}/_data/uhf_{s}_{opt}_noons.txt"
        log_file = f"../../uhf_survey_no_fn_dft_diis/{s}/_logs/uhf_{s}_{opt}.out"
        series.append(get_scf_data(log_file, noon_file, s, opt))


df = pd.concat(series, axis=1,).T
print(df)
df.to_csv("_data/uhf_survey_no_fn_dft_diis.csv")
