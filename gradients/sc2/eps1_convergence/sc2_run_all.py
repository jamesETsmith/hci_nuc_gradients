"""
This script calculates the Sc2 gradients using CASSCF as a function of $\epsilon_1$.
See `README.md` for more details.

Author: James E. T. Smith <james.smith9113@gmail.com>
Date: 2/2/2020
"""
import numpy as np
from pyscf import gto, dft, mcscf, grad
from pyscf.shciscf import shci
from shci_gradients import grad_test_utils as gtu

#
# User Settings
#


do_vHCISCF = 1
do_vHCISCF_AA = 1
do_HCISCF = 1
do_HCISCF_AA = 1

#
# Data generating functions
#


def make_sc2(output=None):
    """
    Set up the Sc2 molecule.
    
    Parameters
    ----------
    output : str, optional
        Output file path, by default None
    
    Returns
    -------
    mc
        `pyscf.mcscf.CASSCF` object
    """

    mol = gto.Mole()
    # In Angstroms
    mol.atom = """Sc 0 0 0; Sc 0 0 4.5;"""
    mol.basis = "ccpvdz"
    mol.unit = "Bohr"
    mol.verbose = 4
    mol.spin = 4
    mol.verbose = 5
    mol.symmetry = True
    if output is not None:
        mol.output = output
    mol.build()

    mf = dft.UKS(mol, xc="bp86")
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-10
    mf.kernel()

    ncas = 18
    nelecas = 6

    mc = mcscf.CASSCF(mf, ncas, nelecas)
    mc.fcisolver.wfnsym = "A2u"
    mc.conv_tol = 1e-12
    mc.conv_tol_grad = 1e-6
    return mc


def vHCISCF(epsilon1):
    mc2 = make_sc2("_logs/vhciscf_{}.out".format(epsilon1))
    mc2.fcisolver = shci.SHCI(mc2.mol)
    mc2.fcisolver.useExtraSymm = True
    mc2.fcisolver.sweep_epsilon = [epsilon1]
    mc2.fcisolver.sweep_iter = [0]
    mc2.mc2step()

    grad2 = mc2.Gradients().kernel()
    np.savetxt("_data/vhciscf_{}.txt".format(epsilon1), grad2)


def vHCISCF_AA(epsilon1):
    mc3 = make_sc2("_logs/vhciscf_aa_{}.out".format(epsilon1))
    mc3.fcisolver = shci.SHCI(mc3.mol)
    mc3.fcisolver.useExtraSymm = True
    mc3.fcisolver.sweep_epsilon = [epsilon1]
    mc3.fcisolver.sweep_iter = [0]
    mc3.mc2step()

    mc3.internal_rotation = True
    mc3.max_cycle_macro = 200
    mc3.mc2step()

    grad3 = mc3.Gradients().kernel()
    np.savetxt("_data/vhciscf_aa_{}.txt".format(epsilon1), grad3)


def HCISCF(epsilon1):
    mc4 = make_sc2("_logs/hciscf_{}.out".format(epsilon1))
    mc4.fcisolver = shci.SHCI(mc4.mol)
    mc4.fcisolver.useExtraSymm = True
    mc4.fcisolver.sweep_epsilon = [epsilon1]
    mc4.fcisolver.sweep_iter = [0]
    mc4.fcisolver.stochastic = False
    mc4.fcisolver.epsilon2 = 1e-10
    mc4.mc2step()

    grad4 = mc4.Gradients().kernel()
    np.savetxt("_data/hciscf_{}.txt".format(epsilon1), grad4)


def HCISCF_AA(epsilon1):
    mc5 = make_sc2("_logs/hciscf_aa_{}.out".format(epsilon1))
    mc5.fcisolver = shci.SHCI(mc5.mol)
    mc5.fcisolver.useExtraSymm = True
    mc5.fcisolver.sweep_epsilon = [epsilon1]
    mc5.fcisolver.sweep_iter = [0]
    mc5.fcisolver.stochastic = False
    mc5.fcisolver.epsilon2 = 1e-10
    mc5.mc2step()

    mc5.internal_rotation = True
    mc5.max_cycle_macro = 200
    mc5.mc2step()

    grad5 = mc5.Gradients().kernel()
    np.savetxt("_data/hciscf_aa_{}.txt".format(epsilon1), grad5)


if __name__ == "__main__":
    epsilon1s = [5e-3, 1e-3, 1e-4, 1e-5, 1e-6]

    if do_vHCISCF:
        for eps1 in epsilon1s:
            vHCISCF(eps1)

    if do_vHCISCF_AA:
        for eps1 in epsilon1s:
            vHCISCF_AA(eps1)

    if do_HCISCF:
        for eps1 in epsilon1s:
            HCISCF(eps1)

    if do_HCISCF_AA:
        for eps1 in epsilon1s:
            HCISCF_AA(eps1)
