import sys
import os
import numpy as np
import pandas as pd
from rmsd import kabsch_rmsd

sys.path.append("..")
from analysis_utils import tabulate_geometry_errors, parse_multi_XYZ, format_exponential

os.makedirs("_data", exist_ok=True)


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
xyz_files = [
    f"../../dz/{d}/40e_40o/_data/{d}_DZ_40e_40o_eps1=7.5e-05.xyz" for d in states
]


# compare_geom_from_file(xyz_files[0], xyz_files[1], "(1,A)/(1,B)")
# compare_geom_from_file(xyz_files[3], xyz_files[4], "(3,A)/(3,B)")

for i, si in enumerate(states):

    df_i, df_f = tabulate_geometry_errors(xyz_files[i])
    df_i["Geometry"] = ["Initial"] * len(df_i.index)
    df_f["Geometry"] = ["Final"] * len(df_f.index)
    df = pd.concat([df_i, df_f])
    df.to_csv(f"_data/{si}.csv")


#
# Compare singlets amongst themselves
#
starting_geoms = [parse_multi_XYZ(f)[1][0] for f in xyz_files if "1_" in f]
final_geoms = [parse_multi_XYZ(f)[1][-1] for f in xyz_files if "1_" in f]
final_triplet_geoms = [parse_multi_XYZ(f)[1][-1] for f in xyz_files if "3_" in f]


df = pd.DataFrame(
    {
        "Geometries": ["A,B", "A,C", "B,C"],
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
