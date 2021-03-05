import numpy as np
import pandas as pd
from noci_utils import noci_stability

np.set_printoptions(linewidth=120)

# Lists of chkfiles
chk_1_A = [
    "../../uhf_survey_dft_newton/1_A/_chk/uhf_1_A_newton.chk",
    "../../uhf_survey_no_fn_dft_diis/1_A/_chk/uhf_1_A_newton.chk",
    "../../uhf_survey_no_fn_dft_adiis/1_A/_chk/uhf_1_A_adiis.chk",
]
chk_1_B = [
    "../../uhf_survey_no_fn_dft_adiis/1_B/_chk/uhf_1_B_newton.chk",
    "../../uhf_survey_dft_no_adiis/1_B/_chk/uhf_1_B_newton.chk",
    "../../uhf_survey_dft_newton/1_B/_chk/uhf_1_B_newton.chk",
]
chk_1_C = [
    "../../uhf_survey_no_fn_dft_adiis/1_C/_chk/uhf_1_C_newton.chk",
    "../../uhf_survey/1_C/_chk/uhf_1_C_newton.chk",
    "../../uhf_survey_dft_no_adiis/1_C/_chk/uhf_1_C_newton.chk",
]
chk_3_A = [
    "../../uhf_survey_no_fn_dft_diis/3_A/_chk/uhf_3_A_newton.chk",
    "../../uhf_survey_dft_no_adiis/3_A/_chk/uhf_3_A_newton.chk",
    "../../uhf_survey_no_fn_dft_adiis/3_A/_chk/uhf_3_A_newton.chk",
]
chk_3_B = [
    "../../uhf_survey_no_fn_dft_adiis/3_B/_chk/uhf_3_B_newton.chk",
    "../../uhf_survey/3_B/_chk/uhf_3_B_newton.chk",
    "../../uhf_survey_no_fn_dft_diis/3_B/_chk/uhf_3_B_newton.chk",
]
chk_3_C = ["../../uhf_survey/3_C/_chk/uhf_3_C_newton.chk"]

# NOCI Calcs
pp = True
RTOL = 1e-8

# noci_stability(chk_1_A, outputbase="1_A", partial_projection=pp, RTOL=RTOL)
# noci_stability(chk_1_B, outputbase="1_B", partial_projection=pp, RTOL=RTOL)
# noci_stability(chk_1_C, outputbase="1_C", partial_projection=pp, RTOL=RTOL)
noci_stability(chk_3_A, outputbase="3_A", partial_projection=pp, RTOL=RTOL)
noci_stability(chk_3_B, outputbase="3_B", partial_projection=pp, RTOL=RTOL)
# noci_stability(chk_3_C, outputbase="3_C", partial_projection=pp, RTOL=RTOL)

# Collect data
noci_data_series = []
for state in ["1_A", "1_B", "1_C", "3_A", "3_B", "3_C"]:
    noci_data_series.append(
        pd.read_csv(
            f"_data/noci_summaries/{state}.csv",
            # header=None,
            index_col=0,
            squeeze=True,
        )
    )

df = pd.DataFrame(noci_data_series)
print(df.to_string())
df.to_csv(f"_data/noci_summaries/noci_summary_PP={pp}_RTOL={RTOL}.csv")
