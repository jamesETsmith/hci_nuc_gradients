"""
Calculate the gradient as function of epsilon_1 for vHCISCF.
"""

import numpy as np

from pyscf import gto, scf, mcscf, grad
from pyscf.shciscf import shci


def run_hciscf(eps1: float):
    mol = gto.Mole()
    # In Angstroms
    mol.atom = """N 0 0 0; N 0 0 1;"""
    mol.basis = "ccpvdz"
    mol.verbose = 4
    mol.symmetry = True
    mol.output = "_logs/vhciscf_aa_{}.out".format(eps1)
    mol.build()

    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-10
    mf.kernel()

    ncas = 8
    nelecas = 10

    mc = mcscf.CASSCF(mf, ncas, nelecas)
    mc.conv_tol = 1e-12
    mc.conv_tol_grad = 1e-6

    mc.fcisolver = shci.SHCI(mc.mol)
    mc.fcisolver.sweep_epsilon = [eps1]
    mc.fcisolver.sweep_iter = [0]
    mc.mc2step()

    # AA
    mc.internal_rotation = True
    mc.max_cycle_macro = 100
    mc.conv_tol = 1e-7
    mc.conv_tol_grad = 1e-3
    if eps1 >= 1e-3:
        mc.conv_tol *= 10
        mc.conv_tol_grad *= 2
        mc.max_cycle_macro = 200
    mc.mc2step()

    grad = mc.Gradients().kernel()
    np.savetxt("_data/vhciscf_aa_{}.txt".format(eps1), grad)


if __name__ == "__main__":
    epsilon1 = [5e-3, 1e-3, 1e-4, 1e-5, 1e-6]

    for eps1 in epsilon1:
        run_hciscf(eps1)
