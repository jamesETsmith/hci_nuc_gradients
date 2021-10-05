import numpy as np
import pandas as pd
from scipy.constants import calorie, physical_constants, N_A

ha_2_kcalpmol = physical_constants["Hartree energy"][0] * N_A / (calorie * 1e3)
print(ha_2_kcalpmol)


df_i = pd.read_csv("_data/40e_40o_initial_extrap_comparison.csv")
df_f = pd.read_csv("_data/40e_40o_final_extrap_comparison.csv")

# S-T gaps (for initial and final)
# (1,A) - (3,C) (Comparing to Cramer)
# (1,A) - (3,B)

st1_i = df_i["Linear Fit Int. (Ha)"].iloc[5] - df_i["Linear Fit Int. (Ha)"].iloc[0]
st2_i = df_i["Linear Fit Int. (Ha)"].iloc[4] - df_i["Linear Fit Int. (Ha)"].iloc[0]

st1_f = df_f["Linear Fit Int. (Ha)"].iloc[5] - df_f["Linear Fit Int. (Ha)"].iloc[0]
st2_f = df_f["Linear Fit Int. (Ha)"].iloc[4] - df_f["Linear Fit Int. (Ha)"].iloc[0]
# print(st1_i, st1_f)
# print(st2_i, st2_f)

df = pd.DataFrame(
    {
        "Geometry": ["Initial", "Final"],
        " (3,C)-(1,A)": np.array([st1_i, st1_f]) * ha_2_kcalpmol,
        " (3,B)-(1,A)": np.array([st2_i, st2_f]) * ha_2_kcalpmol,
    }
)
print(df)
print(df.to_latex(index=False, float_format="%.1f"))
