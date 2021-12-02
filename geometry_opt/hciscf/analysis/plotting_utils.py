import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from cycler import cycler

# Plotting settings
# From: https://personal.sron.nl/~pault/#sec:qualitative
colorblind_palette = np.array(
    [
        "#EE7733",
        "#0077BB",
        "#33BBEE",
        "#EE3377",
        "#CC3311",
        "#009988",
        "#BBBBBB",
    ]
)

# mpl.rcParams["axes.prop_cycle"] = cycler("color", colorblind_palette)
# sns.set_palette(colorblind_palette)


def set_palette(n_colors: int = 7):
    # Dictionary to store some preset subsets of the full palette
    subset = {
        7: np.arange(7),
        5: np.array([1, 0, 2, 3, 5]),
        4: np.array([0, 1, 3, 5]),
        3: np.array([0, 1, 3]),
        2: np.array([0, 1]),
    }

    mpl.rcParams["axes.prop_cycle"] = cycler(
        "color", colorblind_palette[subset[n_colors]]
    )


def set_context(context: str, font_scale: float = 1.0):
    # font_dirs = "/usr/share/fonts/"
    # font_files = mpl.font_manager.findSystemFonts(fontpaths=font_dirs)

    # for font_file in font_files:
    #     print(font_file)
    #     print(mpl.font_manager.fontManager.addfont(font_file))

    # supported_vals = ["paper", "talk"]
    # print(mpl.rcParams["font.serif"])

    # Default Params for "paper"
    params = {
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["DejaVu Serif"],
        "font.sans-serif": ["DejaVu Serif"],
        # https://matplotlib.org/stable/tutorials/text/mathtext.html
        "mathtext.default": "regular",
        "mathtext.fontset": "cm",
    }

    if context == "talk":
        sns.set_context(context, font_scale=font_scale)
        params["text.usetex"] = False
        params["font.family"] = "sans-serif"
        params["font.serif"] = ["DejaVu Sans"]
        params["font.sans-serif"] = ["DejaVu Sans"]
        params["mathtext.fontset"] = "dejavusans"
    elif context == "paper":
        sns.set_context("notebook", font_scale=font_scale)  # "paper" is too small
    elif context not in supported_vals:
        msg = f"{context} value for context isn't supported try one of the following: "
        for sv in supported_vals:
            msg += f"{sv}, "
        raise ValueError(msg)

    plt.rcParams.update(params)


if __name__ == "__main__":
    set_palette(5)
    set_context("paper")
    # mpl.rcParams.update({"font.size": 36})

    plt.figure()
    for i in range(10):
        plt.scatter([i], [i])

    plt.xlabel("$\mathrm{E_{SHCI}}$ (Ha)")
    plt.tight_layout()
    plt.savefig("test.png", dpi=600)
