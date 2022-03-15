import sys
import os
import numpy as np
import pandas as pd
from rmsd import kabsch_rmsd

sys.path.append("..")
from analysis_utils import tabulate_geometry_errors
from analysis_utils import tabulate_geometry_errors_ndarray
from analysis_utils import parse_multi_XYZ
from analysis_utils import format_exponential

os.makedirs("_data", exist_ok=True)

dft_geoms = np.load("../../../dft/analysis/_data/geoms.npy")


def compare_geom_from_file(f1, f2, msg):

    i_g_1 = parse_multi_XYZ(f1)[1][0]
    f_g_1 = parse_multi_XYZ(f1)[1][-1]
    i_g_2 = parse_multi_XYZ(f2)[1][0]
    f_g_2 = parse_multi_XYZ(f2)[1][-1]

    print(f"RMSD in initial {msg} {kabsch_rmsd(i_g_1, i_g_2):.2e}")
    print(f"RMSD in final   {msg} {kabsch_rmsd(f_g_1, f_g_2):.2e}")


#
# Collect Data
#
states = ["1_A", "1_B", "1_C", "3_A", "3_B", "3_C"]
mn15_idx = [6, None, None, 7, None, 8]
m062x_idx = [0, None, None, 1, None, 2]
xyz_files = [
    f"../../dz/{d}/40e_40o/_data/{d}_DZ_40e_40o_eps1=7.5e-05.xyz" for d in states
]


# compare_geom_from_file(xyz_files[0], xyz_files[1], "(1,A)/(1,B)")
# compare_geom_from_file(xyz_files[3], xyz_files[4], "(3,A)/(3,B)")

for i, si in enumerate(states):

    df_i, df_f = tabulate_geometry_errors(xyz_files[i])
    df_i["Method"] = ["M06-L"] * len(df_i.index)
    df_f["Method"] = ["vHCISCF"] * len(df_f.index)

    data = [df_i, df_f]

    if mn15_idx[i] is not None:
        mn15_geom = dft_geoms[mn15_idx[i]]
        df_mn15 = tabulate_geometry_errors_ndarray(mn15_geom)
        df_mn15["Method"] = ["MN15"] * len(df_f.index)
        data.append(df_mn15)

    if m062x_idx[i] is not None:
        m062x_geom = dft_geoms[m062x_idx[i]]
        df_m062x = tabulate_geometry_errors_ndarray(m062x_geom)
        df_m062x["Method"] = ["M06-2X"] * len(df_f.index)
        data.append(df_m062x)

    df = pd.concat(data)

    df.to_csv(f"_data/{si}.csv")


#
# Compare singlets amongst themselves
#
starting_geoms = [parse_multi_XYZ(f)[1][0] for f in xyz_files if "1_" in f]
final_geoms = [parse_multi_XYZ(f)[1][-1] for f in xyz_files if "1_" in f]
final_triplet_geoms = [parse_multi_XYZ(f)[1][-1] for f in xyz_files if "3_" in f]


df = pd.DataFrame(
    {
        "Geometries": [
            "\textbf{A},\textbf{B}",
            "\textbf{A},\textbf{C}",
            "\textbf{B},\textbf{C}",
        ],
        "Initial": [
            kabsch_rmsd(starting_geoms[0], starting_geoms[1]),
            kabsch_rmsd(starting_geoms[0], starting_geoms[2]),
            kabsch_rmsd(starting_geoms[2], starting_geoms[1]),
        ],
        "Optimized Singlet": [
            kabsch_rmsd(final_geoms[0], final_geoms[1]),
            kabsch_rmsd(final_geoms[0], final_geoms[2]),
            kabsch_rmsd(final_geoms[2], final_geoms[1]),
        ],
        "Optimized Triplet": [
            kabsch_rmsd(final_triplet_geoms[0], final_triplet_geoms[1]),
            kabsch_rmsd(final_triplet_geoms[0], final_triplet_geoms[2]),
            kabsch_rmsd(final_triplet_geoms[2], final_triplet_geoms[1]),
        ],
    }
)
print(df)
df.to_latex(
    "_data/rmsd_geom_comparison.tex",
    index=False,
    float_format=format_exponential,
    escape=False,
)

#
# The Effect of optimization
#
df = pd.DataFrame(
    {
        "Geometries": [r"\textbf{A}", r"\textbf{B}", r"\textbf{C}"],
        "Singlet RMSD": [
            kabsch_rmsd(final_geoms[0], starting_geoms[0]),
            kabsch_rmsd(final_geoms[1], starting_geoms[1]),
            kabsch_rmsd(final_geoms[2], starting_geoms[2]),
        ],
        "Triplet RMSD": [
            kabsch_rmsd(final_triplet_geoms[0], starting_geoms[0]),
            kabsch_rmsd(final_triplet_geoms[1], starting_geoms[1]),
            kabsch_rmsd(final_triplet_geoms[2], starting_geoms[2]),
        ],
    }
)
print(df)
df.to_latex(
    "_data/rmsd_effect_of_opt.tex",
    index=False,
    float_format=format_exponential,
    escape=False,
)
