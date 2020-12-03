"""
This script compares the finite difference gradients to the analytical ones using for the CASSCF method.
See `README.md` for more details.
"""
import copy
import numpy as np

from pyscf import gto, scf, mcscf
from pyscf.shciscf import shci


def make_n2(hi, mo0=None, dm0=None):
    mol = gto.Mole()
    # In Angstroms
    mol.atom = [["N", np.array([0, 0, 0])], ["N", np.array([0, 0, 1 + hi])]]
    mol.basis = "ccpvdz"
    # mol.unit = "Angstrom"
    mol.verbose = 4
    mol.symmetry = False
    mol.build()

    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-10

    if dm0 is None:
        mf.kernel()
    else:
        mf.scf(dm0)

    ncas = 8
    nelecas = 10

    mc = mcscf.CASSCF(mf, ncas, nelecas)
    mc.fcisolver = shci.SHCI(mol)
    mc.fcisolver.sweep_iter = [0]
    mc.fcisolver.sweep_epsilon = [5e-3]
    mc.conv_tol = 1e-12
    mc.conv_tol_grad = 1e-6

    if mo0 is not None:
        mo = mcscf.project_init_guess(mc, mo0)
        mc.mc2step(mo)
    else:
        mc.mc2step()

    # return copy.copy(mc)
    return mc, mc.mo_coeff, mc.make_rdm1()


if __name__ == "__main__":
    import sys
    import json

    from pyscf import grad
    from pyscf.lib import param
    from shci_gradients import grad_test_utils as gtu

    # Check command line input
    h = 1e-4  # Angstroms
    if len(sys.argv) == 2:
        h = float(sys.argv[1])
    elif len(sys.argv) > 2:
        print(sys.argv)
        raise AssertionError("Too man command line arguments (see above).")

    # Analytical Gradient
    mc, mo0, dm0 = make_n2(0)
    g = mc.Gradients()
    pyscf_grad = g.kernel()

    # Finite Difference
    hlist = [-2.0 * h, -h, h, 2.0 * h]
    fd_energies = [make_n2(h_i, mo0, dm0)[0].e_tot for h_i in hlist]

    print(fd_energies)

    # Dump to JSON
    data = {}
    data["h"] = h
    data["central_dif"] = gtu.two_point_stencil(fd_energies[1:3], h / param.BOHR)
    data["five_pt"] = gtu.five_point_stencil(fd_energies, h / param.BOHR)
    data["analytical"] = pyscf_grad[1][2]
    with open("_data/h={}.json".format(h), "w") as f:
        json.dump(data, f)

    # Dump to STDOUT
    print("PySCF Analytical Gradient:\n\t", pyscf_grad)
    print("5pt Stencil Gradient:\n\t", data["five_pt"])
    print("central_dif Gradient:\n\t", data["central_dif"])
    print(abs(pyscf_grad[0][0] - data["five_pt"]))
