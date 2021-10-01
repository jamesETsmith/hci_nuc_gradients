import os
import sys
import numpy as np
from pyscf import gto, scf, dft
from pyscf.tools import molden
from pyscf.lib import chkfile
from pyscf.scf.uhf import mulliken_spin_pop
from scipy.linalg import eigh
from functools import reduce

#
# User Settings
#

# Read in command line args
if len(sys.argv) != 4:
    raise ValueError("Wrong number of CLIs. Use \n\tpython uhf.py 1 A adiis")

multiplicity = int(sys.argv[1])
geometry = sys.argv[2].upper()
opt_strategy = sys.argv[3]

# Checks
if geometry not in ["A", "B", "C"]:
    raise ValueError("Geometry must be A, B, or C.")
if opt_strategy not in ["adiis", "diis", "newton"]:
    raise ValueError("Optimization strategy must be adiis, diis, or newton")

spin = multiplicity - 1
outputbase = f"uhf_{multiplicity}_{geometry}_{opt_strategy}"

# Set up data directories
data_dirs = ["_logs", "_chk", "_molden", "_data"]
for dd in data_dirs:
    if not os.path.exists(dd):
        os.mkdir(dd)


# Build Molecule
mol = gto.Mole()
mol.atom = f"../../_geometries/{geometry}.xyz"
mol.basis = "ccpvtz"
mol.verbose = 5
mol.output = f"_logs/{outputbase}.out"
mol.max_memory = 400000  # in MB
mol.spin = spin
mol.build()


# Without restarting from UKS we get wrong spin density
mf0 = dft.UKS(mol, xc="PBE0")
mf0.diis = scf.CDIIS()
mf0.max_cycle = 100
mf0.kernel()

# UHF
mf = scf.UHF(mol)
if opt_strategy == "newton":
    mf = mf.newton()
elif opt_strategy == "adiis":
    mf.diis = scf.ADIIS()
    mf.max_cycle = 100
elif opt_strategy == "diis":
    mf.diis = scf.CDIIS()
    mf.max_cycle = 100
mf.chkfile = f"_chk/{outputbase}.chk"
mf.kernel(mf0.make_rdm1())

# Check for stability
new_mo, _ = mf.stability()

stable = False
if (np.array(new_mo) == np.array(mf.mo_coeff)).all():
    stable = True
else:
    dm = mf.make_rdm1(new_mo)
    mf.kernel(dm0=dm)

    new_mo, _ = mf.stability()
    if (np.array(new_mo) == np.array(mf.mo_coeff)).all():
        stable = True

mf.analyze()
scf.uhf.mulliken_spin_pop(mol, mf.make_rdm1(), s=mf.get_ovlp())

# Dump orbitals
molden.from_scf(mf, f"_molden/{outputbase}.molden")

# Dump objects for making NO
dm_a, dm_b = mf.make_rdm1()
S = mf.get_ovlp()

Dm = dm_a + dm_b
A = reduce(np.dot, (S, Dm, S))
w, v = eigh(A, b=S)
noons = np.flip(w)
NOs = np.flip(v, axis=1)

np.savetxt(f"_data/{outputbase}_noons.txt", noons)
np.save(f"_chk/{outputbase}_NOs.npy", NOs)
molden.from_mo(mol, f"_molden/{outputbase}_NOs.molden", NOs)
