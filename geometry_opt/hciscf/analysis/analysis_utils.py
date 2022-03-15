import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
from scipy.optimize import curve_fit
from rmsd import kabsch_rmsd

#
# Tools for geomeTRIC
#
def parse_multi_XYZ(filename: str) -> list:
    """Parse energies and geometries"""

    with open(filename, "r") as f:
        lines = f.readlines()

    n_atoms = int(lines[0].split()[0])
    n_iter = int((len(lines)) / (n_atoms + 2))
    # print(len(lines), n_iter)
    assert (len(lines)) % (n_atoms + 2) == 0

    geometries = []  #  in Bohr
    energies = []  # in Ha
    atoms = [None] * n_atoms

    for i in range(n_iter):
        idx = (n_atoms + 2) * i + 2
        energies.append(float(lines[idx - 1].split()[-1]))
        coords = []
        for a in range(n_atoms):
            lsplit = lines[idx + a].split()
            # print(lsplit)
            atoms[a] = lsplit[0]
            xyz = np.array([float(xi) for xi in lsplit[1:]])
            coords.append(xyz)

        geometries.append(np.array(coords))

    return energies, geometries, atoms


def get_shci_data(filename: str) -> dict:
    data = {"var": 0.0, "pt": 0.0, "pt uncertainty": 1.0}

    if not os.path.exists(filename):
        print(f"WARNING: {filename} not found, returning 0")
        return data

    with open(filename, "r") as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if "PTEnergy:" in line and len(line.split()) == 4:
            data["pt"] = float(line.split()[1])
            data["pt uncertainty"] = float(line.split()[3])

        if "Root" in line and "Davidson" not in line:
            data["var"] = float(lines[i + 1].split()[1])

    if data["pt"] == 0.0 and data["pt uncertainty"] == 1.0:
        print(f"Getting SHCI data from {filename} FAILED")

    return data


def extrapolate(
    data: pd.DataFrame, species: str, cas: str, initial_or_final: str
) -> dict:

    if initial_or_final != "initial" and initial_or_final != "final":
        raise ValueError("initial_or_final must be either 'initial' or 'final'")

    #
    # Fitting with Error
    #
    fit_x = np.linspace(data["pt correction"].min(), 0)

    def linear(x, m, b):
        return m * x + b

    def quadratic(x, c1, c2, c3):
        return c1 * np.power(x, 2) + c2 * x + c3

    lin_popt, lin_pcov = curve_fit(
        linear,
        data["pt correction"],
        data["pt"],
        # p0=[slope, intercept],
        sigma=data["pt uncertainty"],
        absolute_sigma=True,
    )

    lin_intercept = lin_popt[-1]
    lin_int_uncert = np.sqrt(np.diag(lin_pcov))[-1]

    quad_popt, quad_pcov = curve_fit(
        quadratic,
        data["pt correction"],
        data["pt"],
        # p0=q_coeffs,
        sigma=data["pt uncertainty"],
        absolute_sigma=True,
    )
    quad_intercept = quad_popt[-1]
    quad_int_uncert = np.sqrt(np.diag(quad_pcov))[-1]

    print(
        f"Linear Fit Intercept    {lin_intercept:.6f} \u00B1 {lin_int_uncert:.1e} (Ha)"
    )
    print(
        f"Quadratic Fit Intercept {quad_intercept:.6f} \u00B1 {quad_int_uncert:.1e} (Ha)"
    )
    extrap_data = {
        "Multiplicity": species[0],
        "Geometry": species[2],
        "Initial/Final": initial_or_final,
        "CAS": cas,
        "Linear Fit Int. (Ha)": lin_intercept,
        "Linear Fit Int. Uncertainty (Ha)": lin_int_uncert,
        "Quadratic Fit Int. (Ha)": quad_intercept,
        "Quadratic Fit Int. Uncertainty (Ha)": quad_int_uncert,
        "Linear Slope": lin_popt[0],
        # "fit_x": fit_x,
        # "lin_fit": linear(fit_x, lin_popt[0], lin_popt[1]),
        # "quad_fit": quadratic(fit_x, quad_popt[0], quad_popt[1], quad_popt[2]),
    }

    #
    # Plot
    #
    plt.figure()
    sns.set_style("ticks")
    sns.set_palette("muted")
    sns.set_context("talk")

    # plt.errorbar(
    #     data["pt correction"], data["pt"], yerr=data["pt uncertainty"], label="Data"
    # )
    # Error bars here are smaller than the circles =)
    plt.plot(data["pt correction"], data["pt"], "o", label="Data")
    plt.plot(fit_x, linear(fit_x, lin_popt[0], lin_popt[1]), label="Linear Fit")
    plt.plot(
        fit_x,
        quadratic(fit_x, quad_popt[0], quad_popt[1], quad_popt[2]),
        label="Quad. Fit",
    )
    plt.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))

    plt.xlabel("$E_2$ (Ha)")
    plt.ylabel(r"E$_{SHCI}$ (Ha)")
    plt.title(f"{species} {cas} {initial_or_final}")
    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
    plt.gca().get_yaxis().get_major_formatter().set_scientific(False)
    plt.legend()
    plt.tight_layout()
    os.makedirs("_figures", exist_ok=True)
    plt.savefig(f"_figures/{species}_{cas}_{initial_or_final}.png", dpi=600)
    plt.savefig(f"_figures/{species}_{cas}_{initial_or_final}.pdf")

    return extrap_data


#
# Geometry comparison helpers and Stieber Data
#
stieber_data = [
    # In angstroms
    ["Fe(1)−N(1)", 1.914, 0, 2],
    ["Fe(1)−N(2)", 1.842, 0, 1],
    ["Fe(1)−N(3)", 1.904, 0, 3],
    ["Fe(1)−N(4)", 1.799, 0, 4],
    ["N(4)−N(5)", 1.117, 4, 5],
    ["N(1)−C(2)", 1.344, 1, 6],
    ["N(3)−C(8)", 1.342, 3, 9],
    ["N(2)−C(3)", 1.370, 2, 7],
    ["N(2)−C(7)", 1.372, 2, 8],
    ["C(2)−C(3)", 1.440, 6, 7],
    ["C(7)−C(8)", 1.432, 8, 9],
    # In Degrees
    ["N(1)−Fe(1)−N(2)", 80.82, 1, 0, 2],
    ["N(1)−Fe(1)−N(3)", 160.71, 1, 0, 3],
    ["N(1)−Fe(1)−N(4)", 98.81, 1, 0, 4],
    ["N(2)−Fe(1)−N(3)", 80.64, 2, 0, 3],
    ["N(2)−Fe(1)−N(4)", 169.42, 2, 0, 4],
    ["N(3)−Fe(1)−N(4)", 98.44, 3, 0, 4],
]


def calculate_bl(indices: list, geom: np.ndarray) -> float:
    atom1 = geom[indices[0], :]
    atom2 = geom[indices[1], :]
    return np.linalg.norm(np.array(atom1) - np.array(atom2))


def calculate_bond_angle(indices: list, geom: np.ndarray) -> float:
    atom1 = geom[indices[0], :]
    center = geom[indices[1], :]
    atom2 = geom[indices[2], :]
    v1 = np.array(atom1) - np.array(center)
    v2 = np.array(atom2) - np.array(center)
    alpha = np.arccos(v1.dot(v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    return 180 * alpha / np.pi


def tabulate_geometry_errors(path: str):
    """Return the errors in the initial and final geometry.

    Parameters
    ----------
    path : str
        [description]
    """

    _, geometries, _ = parse_multi_XYZ(path)

    initial_geom = geometries[0]
    final_geom = geometries[-1]
    print(
        f"Difference between first and last geometry {kabsch_rmsd(initial_geom, final_geom)}"
    )

    initial_data = {"Bond": [], "Error": []}
    final_data = {"Bond": [], "Error": []}
    for dp in stieber_data:
        if len(dp) == 4:
            # Calculate error at initial
            bl = calculate_bl(dp[2:], initial_geom)
            initial_data["Bond"].append(dp[0])
            initial_data["Error"].append(abs(bl - dp[1]))

            # Calculate error at final
            bl = calculate_bl(dp[2:], final_geom)
            final_data["Bond"].append(dp[0])
            final_data["Error"].append(abs(bl - dp[1]))

    return pd.DataFrame(initial_data), pd.DataFrame(final_data)


def tabulate_geometry_errors_ndarray(geom: np.ndarray):
    data = {"Bond": [], "Error": []}
    final_data = {"Bond": [], "Error": []}
    for dp in stieber_data:
        if len(dp) == 4:  # Just get bond lengths
            bl = calculate_bl(dp[2:], geom)
            data["Bond"].append(dp[0])
            data["Error"].append(abs(bl - dp[1]))
    return pd.DataFrame(data)


def parse_mcscf_energies(filename: str) -> np.ndarray:
    with open(filename, "r") as f:
        lines = f.readlines()

    mcscf_energies = []

    for i, li in enumerate(lines):
        if "CASSCF E" in li:
            mcscf_energies.append(float(li.split()[10]))

    return np.array(mcscf_energies)


def get_natural_orbital_energy(filename: str) -> float:
    with open(filename, "r") as f:
        lines = f.readlines()

    for i, li in enumerate(lines):
        if "In Natural orbital, CASCI energy" in li:
            energy = float(li.split()[-1])
            return energy


def get_hci_size(filename: str) -> int:

    with open(filename, "r") as f:
        lines = f.readlines()

    for i, li in enumerate(lines):
        # Don't break here in case we want to use this on an MCSCF output
        if "Performing final tight davidson with tol" in li:
            n_dets = int(lines[i - 1].split()[3])
            # print(n_dets)

    return n_dets


#
# Helpers for making tables
#
def format_exponential(x: float) -> str:
    """Format exponential float values to look pretty in LaTeX."""
    if x == 0:
        return "0.0"

    float_str = f"{x:e}"
    base, exponent = float_str.replace("e-0", "e-").replace("e+0", "e+").split("e")
    return r"${0:.1f} \cdot 10^{{{1}}}$".format(float(base), exponent)
