"""
Utilities for Comparing Analytical N2 Gradients
"""

from pyscf import gto, scf, mcscf, grad


def make_n2(output: str = None) -> mcscf.CASSCF:
    """
    Set up the N2 molecule.
    
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
    mol.atom = """N 0 0 0; N 0 0 1;"""
    mol.basis = "ccpvdz"
    mol.verbose = 4
    mol.symmetry = True
    if output is not None:
        mol.output = output
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
    return mc
