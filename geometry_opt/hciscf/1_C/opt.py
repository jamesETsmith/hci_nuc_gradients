"""
This script runs an HCISCF geometry optimization on Fe(PDI) at the singlet 
starting geometry as described in the paper.

The outline for this script is:

1) Run single-point HCISCF calculation and save integrals as FCIDUMP_INITIAL
    so we can use them to calculate vertical excitation energies.

2) Run a geometry optimization.

3) At the final geometry, run a single-point HCISCF calculation, convert the
    MR MOs to NOs and save them as FCIDUMP_FINAL.

4) Save the NOONs in `_data`.

Author: James E T Smith <james.smith9113@gmail.com>
Date: 4/28/2020
"""
import os, sys
import numpy as np

from pyscf import gto, scf, dft, mcscf, grad, lib
from pyscf.lib import chkfile
from pyscf.shciscf import shci
from pyscf.tools import molden

#
# Read in command line args for epsilon1 and active space
#
eps1 = None
ncas = None

if len(sys.argv) != 3:
    raise ValueError(
        "Wrong number of commandline args, use the following syntax\n\t"
        + "python opt.py 7.5e-5 10"
    )
else:
    eps1 = float(sys.argv[1])
    ncas = int(sys.argv[2])
    print(f"Setting eps1 = {eps1}")
    print(f"Setting ncas = {ncas}")

#
# Setting active space and outputbase
#
nelecas = ncas
outputbase = f"1_C_{nelecas}e_{ncas}o_eps1={eps1}"
scratchdir = os.environ.get("TMP")
nprocs = os.environ.get("OMP_NUM_THREADS")  # This is used for MPI in Dice

uhf_chkfile = "../gen_uno/_chk/uhf_1_C.chk"
uno_chkfile = "../gen_uno/_chk/1_C_100e_100o_eps1=0.001.chk"

#
# Read in DFT checkpoint file
#

uno_mo_guess = chkfile.load(uno_chkfile, "mcscf/mo_coeff")

#
# MCSCF Setup
#
mol = chkfile.load_mol(uno_chkfile)
mol.max_memory = 20000
mol.verbose = 6

# MF
mf = scf.UHF(mol)
mf.__dict__.update(chkfile.load(uhf_chkfile, "scf"))

# Building CASSCF Object
mc = shci.SHCISCF(mf, ncas, nelecas).density_fit()
mc.chkfile = f"_chk/{outputbase}.chk"
mc.fcisolver.sweep_iter = [0, 9, 15]
mc.fcisolver.sweep_epsilon = [eps1] * 3
mc.fcisolver.mpiprefix = f"mpirun -np {nprocs}"
mc.fcisolver.scratchDirectory = f"{scratchdir}/fepdi/singlet/{outputbase}"
mc.fcisolver.integralFile = "FCIDUMP_INITIAL"
if not os.path.exists(mc.fcisolver.scratchDirectory):
    os.makedirs(mc.fcisolver.scratchDirectory)

mc.conv_tol = 1e-6
mc.natorb = True
mc.mc1step(uno_mo_guess)
np.savetxt(f"_data/{outputbase}_noons_initial.txt", mc.mo_occ)


#
# Geometry optimization
#


def opt_callback(ld: dict):

    # Grab the MC object inside the gradient scanner
    mc = ld["g_scanner"].base
    i = ld["self"].cycle

    data_dir = os.path.join("_data", f"{outputbase}_geoms")
    if not os.path.isdir(data_dir):
        os.makedirs(data_dir)
    mc.mol.tofile(os.path.join(data_dir, f"{outputbase}_geom_{i}.xyz"))

    # Fake convergence so opt doesn't fail
    if not mc.converged:
        msg = "MC object not converged, overriding and continuing optimization"
        lib.logger.info(ld["mol"], msg)
        mc.converged = True
        return
    else:
        return


# Convergence threshholds for the optimization
scale = 2.0
conv_params = {
    "convergence_energy": scale * 1e-6,  # Eh
    "convergence_grms": scale * 3e-4,  # Eh/Bohr
    "convergence_gmax": scale * 4.5e-4,  # Eh/Bohr
    "convergence_drms": scale * 1.2e-3,  # Angstrom
    "convergence_dmax": scale * 1.8e-3,  # Angstrom
    "prefix": f"_data/{outputbase}",
}

# Run geom. opt.
mc.natorb = False
mc.fcisolver.integralFile = "FCIDUMP"
opt = mc.Gradients().optimizer(solver="geomeTRIC")
opt.callback = opt_callback
mol_eq = opt.kernel(conv_params)

chkfile.save_mol(mol_eq, f"_chk/{outputbase}_mol.chk")
np.save(f"_chk/{outputbase}.npy", mc.mo_coeff)

# Run calc at final geometry to get NOONs and FCIDUMP
mf2 = scf.RHF(mol_eq)
mf2.mo_coeff = mc.mo_coeff

mc = shci.SHCISCF(mf2, ncas, nelecas).density_fit()
mc.chkfile = f"_chk/{outputbase}.chk"
mc.fcisolver.sweep_iter = [0, 9, 15]
mc.fcisolver.sweep_epsilon = [eps1] * 3
mc.fcisolver.mpiprefix = f"mpirun -np {nprocs}"
mc.fcisolver.scratchDirectory = f"{scratchdir}/fepdi/singlet/{outputbase}"
mc.fcisolver.integralFile = "FCIDUMP_FINAL"

mc.conv_tol = 1e-6
mc.natorb = True
mc.mc1step(chkfile.load(mc.chkfile, "mcscf/mo_coeff"))


# Save Final NOs and NOONs
np.save(f"_chk/{outputbase}_NO.npy", mc.mo_coeff)
molden.from_mo(mol_eq, f"_molden/{outputbase}.molden", mc.mo_coeff)
np.savetxt(f"_data/{outputbase}_noons_final.txt", mc.mo_occ)

# 
mc.analyze()