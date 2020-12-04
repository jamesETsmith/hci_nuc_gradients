"""
This script calculates the Sc2 gradients using CASSCF.
See `README.md` for more details.

Author: James E. T. Smith <james.smith9113@gmail.com>
Date: 2/7/2020
"""
import numpy as np
from pyscf import gto, dft, mcscf, grad
from pyscf.shciscf import shci
from shci_gradients import grad_test_utils as gtu


def make_sc2(output: str = None) -> mcscf.CASSCF:
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


#
# User Settings
#


CASSCF = 1
vHCISCF = 1
vHCISCF_AA = 1
HCISCF = 1
HCISCF_AA = 1

epsilon1 = 5e-3

#
# CASSCF
#
if CASSCF:
    mc = make_sc2("_logs/casscf.out")
    mc.mc2step()

    grad = mc.Gradients().kernel()
    np.save("_data/grad_casscf", grad)
    del mc

#
# vHCISCF
#
if vHCISCF:
    mc2 = make_sc2("_logs/vhciscf.out")
    mc2.fcisolver = shci.SHCI(mc2.mol)
    mc2.fcisolver.useExtraSymm = True
    mc2.fcisolver.sweep_epsilon = [epsilon1]
    mc2.fcisolver.sweep_iter = [0]
    mc2.mc2step()

    grad2 = mc2.Gradients().kernel()
    np.save("_data/grad_vhciscf", grad2)
    del mc2

#
# vHCISCF + AA
#
if vHCISCF_AA:
    mc3 = make_sc2("_logs/vhciscf_aa.out")
    mc3.fcisolver = shci.SHCI(mc3.mol)
    mc3.fcisolver.useExtraSymm = True
    mc3.fcisolver.sweep_epsilon = [epsilon1]
    mc3.fcisolver.sweep_iter = [0]
    mc3.conv_tol = 1e-8
    mc3.conv_tol_grad = 2e-4
    mc3.mc2step()

    mc3.internal_rotation = True
    mc3.max_cycle_macro = 200
    mc3.conv_tol = 1e-6
    mc3.conv_tol_grad = 3e-6
    mc3.mc2step()

    grad3 = mc3.Gradients().kernel()
    np.save("_data/grad_vhciscf_aa", grad3)

    del mc3


#
# HCISCF
#
if HCISCF:
    mc4 = make_sc2("_logs/hciscf.out")
    mc4.fcisolver = shci.SHCI(mc4.mol)
    mc4.fcisolver.useExtraSymm = True
    mc4.fcisolver.sweep_epsilon = [epsilon1]
    mc4.fcisolver.sweep_iter = [0]
    mc4.fcisolver.stochastic = False
    mc4.fcisolver.epsilon2 = 1e-10
    mc4.mc2step()

    grad4 = mc4.Gradients().kernel()
    np.save("_data/grad_hciscf", grad4)

    del mc4
#
# HCISCF + AA
#
if HCISCF_AA:
    mc5 = make_sc2("_logs/hciscf_aa.out")
    mc5.fcisolver = shci.SHCI(mc5.mol)
    mc5.fcisolver.useExtraSymm = True
    mc5.fcisolver.sweep_epsilon = [epsilon1]
    mc5.fcisolver.sweep_iter = [0]
    mc5.fcisolver.stochastic = False
    mc5.fcisolver.epsilon2 = 1e-10
    mc5.conv_tol = 1e-8
    mc5.conv_tol_grad = 2e-4
    mc5.mc2step()

    mc5.internal_rotation = True
    mc5.max_cycle_macro = 200
    mc5.conv_tol = 1e-6
    mc5.conv_tol_grad = 3e-6
    mc5.mc2step()

    grad5 = mc5.Gradients().kernel()
    np.save("_data/grad_hciscf_aa", grad5)

    del mc5
