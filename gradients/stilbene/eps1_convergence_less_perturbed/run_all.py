"""
This script calculates the stilbene gradients using CASSCF.
See `README.md` for more details.

Author: James E. T. Smith <james.smith9113@gmail.com>
Date: 3/30/2022
"""
import os
import numpy as np
from pyscf import gto, scf, dft, mp, mcscf, grad
from pyscf.shciscf import shci



def run_stilbene(eps1:float, hci_type:str, aa:bool):
    """
    Set up the stilbene molecule.

    Parameters
    ----------
    output : str, optional
        Output file path, by default None

    Returns
    -------
    mc
        `pyscf.mcscf.CASSCF` object
    """


    if aa:
        outputbase = f"{hci_type.lower()}_aa_{eps1:.1e}"
    else:
        outputbase = f"{hci_type.lower()}_{eps1:.1e}"

    print(outputbase)

    mol = gto.Mole()
    # In Angstroms
    mol.atom = "stilbene.xyz"
    mol.basis = "ccpvdz"
    mol.spin = 0
    mol.verbose = 6
    mol.symmetry = True
    mol.max_memory = 100000
    mol.output = f"_logs/{outputbase}.out"
    mol.build()

    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-10
    mf.kernel()

    mymp = mp.MP2(mf).run()
    noons, natorbs = mcscf.addons.make_natural_orbitals(mymp)
    print(noons)
    print(noons[(noons > 0.02) & (noons < 1.98)])

    ncas = 14
    nelecas = 14

    mf.mo_coeff = natorbs
    mc = mcscf.CASSCF(mf, ncas, nelecas)
    mc.conv_tol = 1e-10
    mc.conv_tol_grad = 1e-5

    if hci_type != "CASSCF":
        mc.fcisolver = shci.SHCI(mc.mol)
        mc.fcisolver.sweep_epsilon = [eps1]*3
        mc.fcisolver.sweep_iter = [0, 6, 12]
        mc.fcisolver.scratchDirectory = "/tmp"
        mc.fcisolver.runtimeDir = "/tmp"
        mc.fcisolver.mpiprefix = "mpirun -np 64"

    if hci_type == "vHCISCF" and not aa:
        mc.mc2step()

    elif hci_type == "vHCISCF" and aa:
        mc.conv_tol = 1e-8
        mc.conv_tol_grad = 2e-4
        mc.mc2step()

        mc.internal_rotation = True
        mc.max_cycle_macro = 200
        mc.conv_tol = 1e-10
        mc.conv_tol_grad = 1e-4
        mc.mc2step()

    elif hci_type == "HCISCF" and not aa:
        mc.fcisolver.stochastic = False
        mc.fcisolver.epsilon2 = 1e-12
        mc.mc2step()

    elif hci_type == "HCISCF" and aa:
        mc.fcisolver.stochastic = False
        mc.fcisolver.epsilon2 = 1e-12
        mc.conv_tol = 1e-8
        mc.conv_tol_grad = 2e-4
        mc.mc2step()

        mc.internal_rotation = True
        mc.max_cycle_macro = 200
        mc.conv_tol = 1e-10
        mc.conv_tol_grad = 1e-4
        mc.mc2step()
        
    elif hci_type == "CASSCF":
        mc.mc2step()

    grad = mc.Gradients().kernel()
    np.save(f"_data/{outputbase}", grad)



if __name__ == "__main__":
    import argparse

    # fmt: off
    parser = argparse.ArgumentParser(description="""This script runs an vHCISCF gradient calculations.""")
    parser.add_argument("--eps1", help="Epsilon 1 for the HCI calculation", type=float, required=True)
    parser.add_argument("--hci_type", help="The type of HCI to use, i.e. vHCISCF, HCISCF", type=str, choices=["vHCISCF", "HCISCF", "CASSCF"], required=True)
    parser.add_argument("--aa", help="Whether to use active-active rotations.", default=False, action="store_true")
    args = parser.parse_args()
    # fmt: on

    os.makedirs("_logs", exist_ok=True)
    os.makedirs("_data", exist_ok=True)

    run_stilbene(args.eps1, args.hci_type, args.aa)
