import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.optimize import curve_fit


def get_shci_data(filename: str) -> dict:
    data = {"var": 0, "pt": 0, "pt uncertainty": 1}

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
        "fit_x": fit_x,
        "lin_fit": linear(fit_x, lin_popt[0], lin_popt[1]),
        "quad_fit": quadratic(fit_x, quad_popt[0], quad_popt[1], quad_popt[2]),
    }
    

    # #
    # # Plot
    # #
    # plt.figure()
    # sns.set_style("ticks")
    # sns.set_palette("muted")

    # # plt.errorbar(
    # #     data["pt correction"], data["pt"], yerr=data["pt uncertainty"], label="Data"
    # # )
    # # Error bars here are smaller than the circles =)
    # plt.plot(data["pt correction"], data["pt"], "o", label="Data")
    # plt.plot(fit_x, linear(fit_x, lin_popt[0], lin_popt[1]), label="Linear Fit")
    # plt.plot(
    #     fit_x,
    #     quadratic(fit_x, quad_popt[0], quad_popt[1], quad_popt[2]),
    #     label="Quad. Fit",
    # )
    # plt.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))

    # plt.xlabel("$E_2$ (Ha)")
    # plt.ylabel(r"E$_{SHCI}$ (Ha)")
    # plt.title(f"{species} {cas} {initial_or_final}")
    # plt.legend()
    # plt.tight_layout()
    # os.makedirs("_figures", exist_ok=True)
    # plt.savefig(f"_figures/{species}_{cas}_{initial_or_final}.png", dpi=600)
    # plt.savefig(f"_figures/{species}_{cas}_{initial_or_final}.pdf")

    return extrap_data
