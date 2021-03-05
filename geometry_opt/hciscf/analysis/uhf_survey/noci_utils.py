import os
import scipy
import numpy as np
import pandas as pd
from pyscf import gto, scf, cc
from pyscf.lib import chkfile


def get_orthoAO(S, LINDEP_CUTOFF=1e-10):
    sdiag, Us = np.linalg.eigh(S)
    print(
        f"Overlap: {sdiag[sdiag < LINDEP_CUTOFF].size}/{sdiag.size} eigenvalues are 0"
    )
    # _, s, _ = np.linalg.svd(S)
    # print(f"Eigenvalues of S {sdiag}")
    # print(f"Singular Values of S {s}")

    X = Us[:, sdiag > LINDEP_CUTOFF] / np.sqrt(sdiag[sdiag > LINDEP_CUTOFF])
    return X


# Green's functions
def gab(A, B, S, RTOL=1e-6):  # S is metric
    ASB = (A.conj().T).dot(S).dot(B)

    # Numerical stability checks
    cond = np.linalg.cond(ASB)
    message = "\tCondition number of (A^dag S B) "
    if cond > 1e2:
        print(message + "\033[91m" + f" {cond:.1e} " + "\033[00m")
    else:
        print(message + "\033[92m" + f" {cond:.1e} " + "\033[00m")

    _, s, _ = np.linalg.svd(ASB)
    print("\tSmallest singular values", s[-2:])
    print(
        f"\tKeeping {s[s>s[0]*RTOL].size} of {s.size} singular vectors (RTOL={RTOL:.1e})"
    )

    if s[-1] < 1e-10:
        raise AssertionError(
            "Smallest singular value of (A^dag S B) in Gab is too small"
        )

    # Actual calc
    inv_O = scipy.linalg.pinv(ASB, rcond=RTOL)
    GAB = B.dot(inv_O.dot(A.conj().T))
    return GAB


def load_from_chk(chk):
    mol = chkfile.load_mol(chk)
    mf = scf.UHF(mol)
    mf.__dict__.update(chkfile.load(chk, "scf"))
    return mol, mf


def noci(chk1, chk2, nstates=1):
    print("#" * 80)
    print("Running NOCI test to determine stability of UHF solutions.\n")
    print(f"Loading 1st UHF from chk={chk1}")
    print(f"Loading 2nd UHF from chk={chk2}")
    # Load data from chkfiles
    mol1, mf1 = load_from_chk(chk1)
    mol2, mf2 = load_from_chk(chk2)
    print(f"Energy from 1st UHF solution: {mf1.e_tot}")
    print(f"Energy from 2nd UHF solution: {mf2.e_tot}")
    # Unpack the two mean-field objects
    mol = mf1.mol
    mo_1, mo_2 = (mf1.mo_coeff, mf2.mo_coeff)
    e_1, e_2 = (mf1.e_tot, mf2.e_tot)
    Sao = mol.intor("int1e_ovlp")
    nocca, noccb = mol.nelec
    print(f"# alpha electrons = {nocca}")
    print(f"# beta  electrons = {noccb}")
    # MO Coefficients
    C1a = mo_1[0]
    C1b = mo_1[1]
    C2a = mo_2[0]
    C2b = mo_2[1]
    C3a = mo_1[1]
    C3b = mo_1[0]
    C4a = mo_2[1]
    C4b = mo_2[0]
    psias = [C1a, C2a, C3a, C4a]
    psibs = [C1b, C2b, C3b, C4b]
    ehfs = [e_1, e_2, e_1, e_2]

    # Construct overlap matrix between solutions
    print(
        "Constructing the NOCI overlap matrix and Hamiltonian for nstates = {}".format(
            nstates
        )
    )
    S = np.zeros((nstates, nstates))
    for i in range(nstates):
        S[i, i] = 1.0  # assume normalized psis
        for j in range(i + 1, nstates):
            S12a = psias[i][:, :nocca].T.dot(Sao).dot(psias[j][:, :nocca])
            S12b = psibs[i][:, :noccb].T.dot(Sao).dot(psibs[j][:, :noccb])
            S[i, j] = np.linalg.det(S12a) * np.linalg.det(S12b)
            S[j, i] = S[i, j]

    # Make Hamiltonian in CI space
    enuc = mf1.energy_nuc()
    H = np.zeros((nstates, nstates))
    for i in range(nstates):
        H[i, i] = ehfs[i]  # assume normalized psis
        for j in range(i + 1, nstates):
            # Off diagonal elements (Green's functions)
            Ga = gab(
                psias[i][:, :nocca], psias[j][:, :nocca], Sao
            )  # alpha spin density matrix (mixed)
            Gb = gab(
                psibs[i][:, :noccb], psibs[j][:, :noccb], Sao
            )  # beta spin density matrix (mixed)
            H12 = (mf1.energy_elec([Ga, Gb])[0] + enuc) * S[i, j]
            H[i, j] = H12
            H[j, i] = H[i, j]
    print()
    print("H =\n{}".format(H))
    print("S =\n{}".format(S))
    X = get_orthoAO(S, 1e-6)
    Ht = X.T.dot(H).dot(X)
    e, v = np.linalg.eigh(Ht)
    idx = e.argsort()
    e = e[idx]
    v = v[:, idx]
    coeffs = X.dot(v[:, 0])  # coeffs in the original basis
    print("NOCI eigenvalues", e)
    print("NOCI ground state coefficients = {}".format(coeffs))
    print("#" * 80 + "\n")

    # 1RDM computation
    nbsf = C1a.shape[0]
    num = np.zeros((nbsf, nbsf))
    denom = 0.0
    for i in range(nstates):
        for j in range(nstates):
            Ga = gab(
                psias[i][:, :nocca], psias[j][:, :nocca], Sao
            )  # alpha spin density matrix (mixed)
            Gb = gab(
                psibs[i][:, :noccb], psibs[j][:, :noccb], Sao
            )  # beta spin density matrix (mixed)
            const = coeffs[i] * coeffs[j] * S[i, j]
            num += const * (Ga + Gb)
            denom += const
    Pao = num / denom
    Xao = get_orthoAO(Sao, 1e-3)
    Xaoinv = np.linalg.inv(Xao)
    Poao = Xaoinv.dot(Pao).dot(Xaoinv.T)
    e, v = np.linalg.eigh(Poao)
    idx = e.argsort()[::-1]
    e = e[idx]
    v = v[:, idx]
    C_natural = Xao.dot(v)
    print("natural orbital occupation numbers : {}".format(e))
    print("nelec = {}".format(np.sum(e)))
    print("nelec = {}".format(Pao.dot(Sao).trace()))

    # <S^2> computation
    Ms = (nocca - noccb) / 2.0
    S2exact = Ms * (Ms + 1.0)
    num = 0.0
    denom = 0.0
    for i in range(nstates):
        for j in range(nstates):
            Ga = gab(
                psias[i][:, :nocca], psias[j][:, :nocca], Sao
            )  # alpha spin density matrix (mixed)
            Gb = gab(
                psibs[i][:, :noccb], psibs[j][:, :noccb], Sao
            )  # beta spin density matrix (mixed)
            const = coeffs[i] * coeffs[j] * S[i, j]
            num += const * (noccb + S2exact - (Gb.dot(Sao).dot(Ga).dot(Sao)).trace())
            denom += const
    S2 = num / denom
    print("<S^2>_NOCI = {}".format(S2))
    return C_natural


def noci_stability(
    chks: str, partial_projection: bool = True, outputbase: str = "", RTOL=1e-3
):
    # Data to return (and optionally dump to csv if outputbase != "")
    data = {
        "State": outputbase,
        "Partial Projection": partial_projection,
        "Gab_RTOL": RTOL,
    }

    # Start the output
    print("#" * 80)
    print("Running NOCI test to determine stability of UHF solutions.\n")
    for chk in chks:
        print(f"Loading UHF from chk={chk}")

    # Load data from chkfiles
    mfs = [load_from_chk(chk)[1] for chk in chks]

    for i, mf in enumerate(mfs):
        print(f"Energy from UHF solution ({i}): {mf.e_tot}")

    # Unpack the common mean-field quantities
    mol = mfs[0].mol
    Sao = mol.intor("int1e_ovlp")
    nmo = Sao.shape[0]
    nocca, noccb = mol.nelec
    enuc = mfs[0].energy_nuc()
    print(f"# alpha electrons = {nocca}")
    print(f"# beta  electrons = {noccb}")

    nstates = len(mfs) * 2 if partial_projection else len(mfs)
    data["nstates"] = nstates
    print(f"Targetting {nstates} states")

    # Unpack MO Coefficients
    psias = [mf.mo_coeff[0] for mf in mfs]
    psibs = [mf.mo_coeff[1] for mf in mfs]
    ehfs = [mf.e_tot for mf in mfs]
    data["Old Energy"] = np.min(ehfs)

    if partial_projection:
        psias += [mf.mo_coeff[1] for mf in mfs]
        psibs += [mf.mo_coeff[0] for mf in mfs]
        ehfs += [mf.e_tot for mf in mfs]

    # Calculate all Green's functions
    print("\nCalculating mixed estimator density matrices")
    Gij_a = np.zeros((nstates, nstates, nmo, nmo))
    Gij_b = np.zeros((nstates, nstates, nmo, nmo))
    for i in range(nstates):
        for j in range(nstates):
            print(f"Calculating Gij between {i} and {j}")
            Gij_a[i, j] = gab(psias[i][:, :nocca], psias[j][:, :nocca], Sao, RTOL=RTOL)
            Gij_b[i, j] = gab(psibs[i][:, :noccb], psibs[j][:, :noccb], Sao, RTOL=RTOL)

    # Construct overlap matrix between solutions
    print("\nConstructing the NOCI overlap matrix and Hamiltonian for nstates")

    S = np.zeros((nstates, nstates))
    for i in range(nstates):
        S[i, i] = 1.0  # assume normalized psis
        for j in range(i + 1, nstates):
            S12a = psias[i][:, :nocca].T.dot(Sao).dot(psias[j][:, :nocca])
            S12b = psibs[i][:, :noccb].T.dot(Sao).dot(psibs[j][:, :noccb])
            S[i, j] = np.linalg.det(S12a) * np.linalg.det(S12b)
            S[j, i] = S[i, j]

    # Make Hamiltonian in CI space
    H = np.zeros((nstates, nstates))
    for i in range(nstates):
        H[i, i] = ehfs[i]  # assume normalized psis
        for j in range(i + 1, nstates):
            # Off diagonal elements (Green's functions)
            # Ga = gab(psias[i][:, :nocca], psias[j][:, :nocca], Sao)
            # Gb = gab(psibs[i][:, :noccb], psibs[j][:, :noccb], Sao)
            Ga, Gb = (Gij_a[i, j], Gij_b[i, j])
            H12 = (mfs[0].energy_elec([Ga, Gb])[0] + enuc) * S[i, j]
            H[i, j] = H12
            H[j, i] = H[i, j]

    print()
    print("H =\n{}".format(H))
    print("S =\n{}".format(S))
    X = get_orthoAO(S, LINDEP_CUTOFF=1e-6)
    print("CHECK DIAG", np.diag(X.T.dot(S.dot(X))))
    Ht = X.T.dot(H).dot(X)
    print(f"Condition number of H {np.linalg.cond(H)}")
    print(f"Condition number of Ht {np.linalg.cond(Ht)}")
    e, v = np.linalg.eigh(Ht)
    idx = e.argsort()
    e = e[idx]
    v = v[:, idx]
    coeffs = X.dot(v[:, 0])  # coeffs in the original basis
    print("NOCI eigenvalues", e)
    print("NOCI ground state coefficients = {}".format(coeffs))
    data["Energy"] = e[0]

    print("\n1-RDM Computation")

    # 1RDM computation
    nbsf = psias[0].shape[0]
    num = np.zeros((nbsf, nbsf))
    denom = 0.0
    for i in range(nstates):
        for j in range(nstates):
            # Ga = gab(psias[i][:, :nocca], psias[j][:, :nocca], Sao)
            # Gb = gab(psibs[i][:, :noccb], psibs[j][:, :noccb], Sao)
            Ga, Gb = (Gij_a[i, j], Gij_b[i, j])
            const = coeffs[i] * coeffs[j] * S[i, j]
            num += const * (Ga + Gb)
            denom += const

    Pao = num / denom
    Xao = get_orthoAO(Sao, LINDEP_CUTOFF=1e-6)
    Xaoinv = np.linalg.inv(Xao)
    Poao = Xaoinv.dot(Pao).dot(Xaoinv.T)
    e, v = np.linalg.eigh(Poao)
    idx = e.argsort()[::-1]
    e = e[idx]
    v = v[:, idx]
    C_natural = Xao.dot(v)

    data["HONO-1"], data["HONO"], data["LUNO"], data["LUNO+1"] = e[117:121]
    if outputbase != "":
        os.makedirs("_data/noci_nos", exist_ok=True)
        np.save(f"_data/noci_nos/{outputbase}.npy", C_natural)
        os.makedirs("_data/noci_noons", exist_ok=True)
        np.savetxt(f"_data/noci_noons/{outputbase}.txt", e)

    print("frontier natural orbital occupation numbers: {}".format(e[113:125]))
    print("nelec = {}".format(np.sum(e)))
    print("nelec = {}".format(Pao.dot(Sao).trace()))

    # <S^2> computation
    print("\nCalculating <S^2>")
    Ms = (nocca - noccb) / 2.0
    S2exact = Ms * (Ms + 1.0)
    num = 0.0
    denom = 0.0
    for i in range(nstates):
        for j in range(nstates):
            # Ga = gab(psias[i][:, :nocca], psias[j][:, :nocca], Sao)
            # Gb = gab(psibs[i][:, :noccb], psibs[j][:, :noccb], Sao)
            Ga, Gb = (Gij_a[i, j], Gij_b[i, j])
            const = coeffs[i] * coeffs[j] * S[i, j]
            num += const * (noccb + S2exact - (Gb.dot(Sao).dot(Ga).dot(Sao)).trace())
            denom += const
    S2 = num / denom
    print("<S^2>_NOCI = {}".format(S2))
    data["S^2"] = S2

    # Dump data
    if outputbase != "":
        os.makedirs("_data/noci_summaries", exist_ok=True)
        pd.Series(data).to_csv(f"_data/noci_summaries/{outputbase}.csv")

    print("#" * 80 + "\n")
    return pd.Series(data)


if __name__ == "__main__":

    # C_no = noci("fe_pdi_uhf_chk/1.chk", "fe_pdi_uhf_chk/2.chk", nstates=4)

    # chk_1_B_diis = "../../uhf_survey_dft_no_adiis/1_B/_chk/uhf_1_B_newton.chk"
    # chk_1_B_newton = "../../uhf_survey_dft_newton/1_B/_chk/uhf_1_B_newton.chk"
    # my_noci(chk_1_B_diis, chk_1_B_newton)
    # exit(0)

    # Check
    # chk_1_B_diis = "../../uhf_survey_dft_no_adiis/1_B/_chk/uhf_1_B_newton.chk"
    # chk_1_B_diis = "../../uhf_survey_dft_no_adiis/1_B/_chk/uhf_1_B_newton.chk"
    # noci(chk_1_B_diis, chk_1_B_diis)
    # exit(0)

    chk_1_B_diis = "../../uhf_survey_dft_no_adiis/1_B/_chk/uhf_1_B_newton.chk"
    chk_1_B_newton = "../../uhf_survey_dft_newton/1_B/_chk/uhf_1_B_newton.chk"
    noci(chk_1_B_diis, chk_1_B_newton, nstates=4)

    noci_stability([chk_1_B_diis, chk_1_B_newton])
    exit(0)

    chk_1_C_adiis = "../../uhf_survey/1_C/_chk/uhf_1_C_newton.chk"
    chk_1_C_diis = "../../uhf_survey_dft_no_adiis/1_C/_chk/uhf_1_C_newton.chk"
    chk_1_C_newton = "../../uhf_survey_dft_newton/1_C/_chk/uhf_1_C_newton.chk"
    noci(chk_1_C_diis, chk_1_C_newton)
