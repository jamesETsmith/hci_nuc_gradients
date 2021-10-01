"""
This script runs an vHCISCF geometry optimization on Fe(PDI) at the singlet 
or triplet starting geometry as described in the paper.

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
import argparse

from pyscf import gto, scf, dft, mcscf, grad, lib
from pyscf.lib import chkfile
from pyscf.shciscf import shci
from pyscf.tools import molden

#
# Read in command line args
#

# fmt: off
parser = argparse.ArgumentParser(description="""This script runs an vHCISCF geometry optimization on Fe(PDI) at the singlet or triplet starting geometry""")
parser.add_argument("--eps1", help="Epsilon 1 for the HCI calculation", type=float, required=True)
parser.add_argument("--ncas", help="The number of active space orbitals, which will be the same number as the electons in the CAS calc.", type=int, required=True)
parser.add_argument("--spin", help="N_alpha - N_beta electrons.", type=int, choices=[0, 2], required=True)
parser.add_argument("--geometry", help="Starting geometry", type=str, choices=["A","B", "C"], required=True)
args = parser.parse_args()
# fmt: on
# Shorthand for some vars
ncas = args.ncas
nelecas = ncas
mult = 2 * (args.spin // 2) + 1
spin = args.spin
eps1 = args.eps1
geom = args.geometry

outputbase = f"{mult}_{geom}_DZ_{nelecas}e_{ncas}o_eps1={eps1}"
print("Output Base:", outputbase)

#
# Setting active space and outputbase
#
scratchdir = os.environ.get("TMP")
nprocs = os.environ.get("OMP_NUM_THREADS")  # This is used for MPI in Dice

# UHF is only used to read in geometries and populate some other mf attributes
# fmt: off
uno_NOs = np.load(f"../gen_uno/_chk/{mult}_{geom}_DZ_100e_100o_eps1=0.0001_vhci_NOs.npy")
uno_mo_occ = np.loadtxt(f"../gen_uno/_data/{mult}_{geom}_DZ_100e_100o_eps1=0.0001_vhci_noons.txt")
# fmt: on


#
# Load the molecule
#
mol = gto.M(
    atom=f"../../../uhf_survey_tz/_geometries/{geom}.xyz",
    basis="ccpvdz",
    spin=spin,
    output=f"_logs/opt_{eps1}.out",
    max_memory=400000,
    verbose=6,
)

# MF
mf = scf.RHF(mol)
mf.mo_coeff = uno_NOs
mf.mo_occ = uno_mo_occ
# mf.__dict__.update(chkfile.load(uhf_chkfile, "scf"))

# Building CASSCF Object
nelecas_by_spin = (nelecas // 2 + spin // 2, nelecas // 2 - spin // 2)
mc = shci.SHCISCF(mf, ncas, nelecas_by_spin).density_fit()
mc.chkfile = f"_chk/{outputbase}.chk"
mc.fcisolver.sweep_iter = [0, 9, 15]
mc.fcisolver.sweep_epsilon = [eps1] * 3
mc.fcisolver.mpiprefix = f"mpirun -np {nprocs}"
mc.fcisolver.scratchDirectory = f"{scratchdir}/fepdi/eps1_response/{outputbase}"
mc.fcisolver.integralFile = f"_molden/FCIDUMP_{outputbase}"
os.makedirs(mc.fcisolver.scratchDirectory, exist_ok=True)

mc.conv_tol = 1e-6
# mc.natorb = True
mc.mc1step(uno_NOs)

# Get gradients and save them
gradients = mc.Gradients().kernel()
np.savetxt(f"_data/{outputbase}_gradients.txt", gradients)

# np.savetxt(f"_data/{outputbase}_noons_initial.txt", mc.mo_occ)
molden.from_mo(mol, f"_molden/{outputbase}_initial.molden", mc.mo_coeff)
mc.analyze()
