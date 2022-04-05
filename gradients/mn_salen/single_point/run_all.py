"""
This script calculates the Mn(salen) gradients using CASSCF.
See `README.md` for more details.

Author: James E. T. Smith <james.smith9113@gmail.com>
Date: 3/30/2022
"""
import os
import numpy as np
from pyscf import gto, scf, dft, mp, mcscf, grad
from pyscf.shciscf import shci


def make_mn_salen(output: str = None) -> mcscf.CASSCF:
    """
    Set up the mn_salen molecule.

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
    mol.atom = "../mn_salen.xyz"
    mol.basis = "ccpvdz"
    mol.spin = 0
    mol.verbose = 6
    mol.symmetry = True
    mol.max_memory = 100000
    if output is not None:
        mol.output = output
    mol.build()

    mf = scf.UHF(mol)
    mf.DIIS = scf.ADIIS
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-6
    mf.max_cycle = 20
    mf.kernel()
    dm1 = mf.make_rdm1()

    mf.DIIS = scf.DIIS
    mf.max_cycle = 50
    mf.kernel(dm0=dm1)

    mo, _, stable, _ = mf.stability(return_status=True)
    for _ in range(5):
        if stable:
            break
        dm1 = mf.make_rdm1(mo, mf.mo_occ)
        mf.kernel(dm0=dm1)
        mo, _, stable, _ = mf.stability(return_status=True)
        mf.mo_coeff = mo



    mymp = mp.MP2(mf).run()
    noons, natorbs = mcscf.addons.make_natural_orbitals(mymp)

    print(noons)
    cas = noons[(noons > 0.02) & (noons < 1.98)]
    print(cas)
    print("NCAS", cas.size)
    print("NELECAS", cas.sum())

    ncas = 22
    nelecas = 28

    mf = scf.addons.convert_to_rhf(mf)
    mf.mo_coeff = natorbs

    mc = mcscf.CASSCF(mf, ncas, (nelecas//2, nelecas//2))
    mc.conv_tol = 1e-10
    mc.conv_tol_grad = 1e-5
    mc.max_cyle_macro = 80
    return mc


#
# User Settings
#
os.makedirs("_logs", exist_ok=True)
os.makedirs("_data", exist_ok=True)

CASSCF = 0
vHCISCF = 1
vHCISCF_AA = 0
HCISCF = 0
HCISCF_AA = 0

epsilon1 = 5e-4

#
# CASSCF
#
if CASSCF:
    mc = make_mn_salen("_logs/casscf.out")
    mc.mc1step()

    grad = mc.Gradients().kernel()
    np.save("_data/grad_casscf", grad)
    del mc

#
# vHCISCF
#
if vHCISCF:
    mc2 = make_mn_salen(f"_logs/vhciscf_{epsilon1:.1e}.out")
    mc2.fcisolver = shci.SHCI(mc2.mol)
    mc2.fcisolver.sweep_epsilon = [epsilon1] *3
    mc2.fcisolver.sweep_iter = [0, 6, 12]
    mc2.fcisolver.scratchDirectory = "/tmp"
    mc2.fcisolver.mpiprefix = "mpirun -np 40"
    mc2.mc1step()

    grad2 = mc2.Gradients().kernel()
    np.save(f"_data/grad_vhciscf_{epsilon1:.1e}", grad2)
    del mc2

#
# vHCISCF + AA
#
if vHCISCF_AA:
    mc3 = make_mn_salen(f"_logs/vhciscf_aa_{epsilon1:.1e}.out")
    mc3.fcisolver = shci.SHCI(mc3.mol)
    mc3.fcisolver.sweep_epsilon = [epsilon1] *3
    mc3.fcisolver.sweep_iter = [0, 6, 12]
    mc3.fcisolver.scratchDirectory = "/tmp"
    mc3.fcisolver.mpiprefix = "mpirun -np 40"
    mc3.conv_tol = 1e-8
    mc3.conv_tol_grad = 2e-4
    mc3.mc1step()

    mc3.internal_rotation = True
    mc3.max_cycle_macro = 200
    mc3.conv_tol = 1e-10
    mc3.conv_tol_grad = 1e-4
    mc3.mc1step()

    grad3 = mc3.Gradients().kernel()
    np.save(f"_data/grad_vhciscf_aa_{epsilon1:.1e}", grad3)

    del mc3


#
# HCISCF
#
if HCISCF:
    mc4 = make_mn_salen(f"_logs/hciscf_{epsilon1:.1e}.out")
    mc4.fcisolver = shci.SHCI(mc4.mol)
    mc4.fcisolver.sweep_epsilon = [epsilon1] *3
    mc4.fcisolver.sweep_iter = [0, 6, 12]
    mc4.fcisolver.scratchDirectory = "/tmp"
    mc4.fcisolver.mpiprefix = "mpirun -np 40"
    mc4.fcisolver.stochastic = False
    mc4.fcisolver.epsilon2 = 1e-12
    mc4.mc1step()

    grad4 = mc4.Gradients().kernel()
    np.save(f"_data/grad_hciscf_{epsilon1:.1e}", grad4)

    del mc4
#
# HCISCF + AA
#
if HCISCF_AA:
    mc5 = make_mn_salen(f"_logs/hciscf_aa_{epsilon1:.1e}.out")
    mc5.fcisolver = shci.SHCI(mc5.mol)
    mc5.fcisolver.sweep_epsilon = [epsilon1] *3
    mc5.fcisolver.sweep_iter = [0, 6, 12]
    mc5.fcisolver.scratchDirectory = "/tmp"
    mc5.fcisolver.mpiprefix = "mpirun -np 40"
    mc5.fcisolver.stochastic = False
    mc5.fcisolver.epsilon2 = 1e-12
    mc5.conv_tol = 1e-8
    mc5.conv_tol_grad = 2e-4
    mc5.mc1step()

    mc5.internal_rotation = True
    mc5.max_cycle_macro = 200
    mc5.conv_tol = 1e-10
    mc5.conv_tol_grad = 1e-4
    mc5.mc1step()

    grad5 = mc5.Gradients().kernel()
    np.save(f"_data/grad_hciscf_aa_{epsilon1:.1e}", grad5)

    del mc5
